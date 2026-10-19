#!/usr/bin/env python3
"""
bench_fit.py -- analysis of bench measurements on a single compass needle
(protocol: docs/BANCADA.md). Extracts the parameters of

    I theta'' = -m B sin(theta - theta0) - b theta' - tau_c sign(theta')

(stuck while |m B sin(theta - theta0)| <= tau_s), i.e. the moment of inertia I,
the magnetic moment m, the viscous damping b, the dynamic (Coulomb) friction
torque tau_c and the static friction torque tau_s. SI units throughout; angles
in radians unless a column name says otherwise.

Subcommands
  geometry  I (and an effective thickness) from the weighed mass and dimensions
  addmass   I from the frequency shift caused by a known added inertia dI
  freq      m/I and the residual field from omega^2 versus coil field
  dipole    m from the on-axis field of the needle, B = (mu0/4pi) 2m / r^3
  decay     b and tau_c from a free decay theta(t): peak-envelope fit, then a
            fit of the full equation of motion to the whole record
  deadband  tau_s from the scatter of the final rest angles of repeated releases
  simulate  synthetic free-decay record (for tests and for planning a run)

Input files are CSV with a header row. Column names are given per subcommand
(--help). Video trackers such as Tracker export t and an angle; use --deg if
the angle is in degrees.
"""
from __future__ import annotations

import argparse
import json
import math
import sys

import numpy as np

MU0_4PI = 1e-7


# ---------------------------------------------------------------------------
# Basic relations
# ---------------------------------------------------------------------------

def inertia_rhombus(mass: float, length: float, width: float,
                    pivot_mass: float = 0.0, pivot_radius: float = 0.0) -> float:
    """Rhombus plate about its centre, I = M (L^2 + W^2) / 24, plus a pivot cylinder."""
    return mass * (length**2 + width**2) / 24.0 + 0.5 * pivot_mass * pivot_radius**2


def inertia_added_mass(f0: float, f1: float, dI: float) -> float:
    """I from the frequencies without (f0) and with (f1) an added inertia dI, same field."""
    if not f1 < f0:
        raise ValueError("the frequency with the added inertia must be lower")
    return dI * f1**2 / (f0**2 - f1**2)


def fit_freq_vs_field(B: np.ndarray, f: np.ndarray):
    """omega^2 = (m/I) (B + B_res): returns m/I, B_res and the fit's relative rms."""
    w2 = (2 * math.pi * np.asarray(f)) ** 2
    A = np.vstack([B, np.ones_like(B)]).T
    (k, c), *_ = np.linalg.lstsq(A, w2, rcond=None)
    rel = float(np.sqrt(np.mean((A @ [k, c] - w2) ** 2)) / np.mean(w2))
    return float(k), float(c / k), rel


def fit_dipole(r: np.ndarray, B: np.ndarray):
    """On-axis dipole field B = 2 mu0/4pi m / r^3 + B_bg: returns m and B_bg."""
    x = 2 * MU0_4PI / np.asarray(r) ** 3
    A = np.vstack([x, np.ones_like(x)]).T
    (m, bg), *_ = np.linalg.lstsq(A, np.asarray(B), rcond=None)
    return float(m), float(bg)


def tau_static_from_deadband(theta_final: np.ndarray, m: float, B: float):
    """Rest angles fill the dead band |theta - theta0| <= arcsin(tau_s / (m B)).

    theta0 is estimated as the centre of the observed range and the half-range
    as the band edge; with n releases the half-range underestimates the edge by
    about 1/(n+1) of the full width, which is corrected for.
    """
    th = np.asarray(theta_final, dtype=float)
    n = th.size
    lo, hi = th.min(), th.max()
    half = 0.5 * (hi - lo) * (n + 1) / max(n - 1, 1)
    return m * B * math.sin(min(half, math.pi / 2)), 0.5 * (hi + lo), half


# ---------------------------------------------------------------------------
# Free decay
# ---------------------------------------------------------------------------

def smooth(t: np.ndarray, th: np.ndarray, period: float) -> np.ndarray:
    """Savitzky-Golay smoothing (cubic) over ~1/4 of a period: tracker noise would
    otherwise add false extrema. The peak bias of this window is < 1e-3 relative."""
    from scipy.signal import savgol_filter
    n = int(round(period / 4.0 / (t[1] - t[0])))
    n = n + 1 if n % 2 == 0 else n
    return savgol_filter(th, n, 3) if n >= 5 else th


def find_extrema(t: np.ndarray, th: np.ndarray, theta0: float | None = None,
                 min_sep: float = 0.0):
    """Alternating extrema of theta - theta0 (parabolic refinement). Returns (t_k, A_k >= 0).

    Pass smoothed data (see smooth()) for noisy records; extrema closer than
    min_sep (use ~T/4) to the previous one are treated as noise.
    """
    if theta0 is None:
        theta0 = float(np.median(th[len(th) // 2:]))
    y = th - theta0
    dy = np.diff(y)
    idx = np.where(np.sign(dy[:-1]) * np.sign(dy[1:]) < 0)[0] + 1
    tk, ak, last_sign = [], [], 0
    for i in idx:
        if i < 1 or i >= len(y) - 1:
            continue
        y0, y1, y2 = y[i - 1], y[i], y[i + 1]
        den = y0 - 2 * y1 + y2
        s = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        s = max(-1.0, min(1.0, s))
        ypk = y1 - 0.25 * (y0 - y2) * s
        sign = 1 if ypk > 0 else -1
        tpk = t[i] + s * (t[i + 1] - t[i])
        if tk and tpk - tk[-1] < min_sep:
            continue
        if sign == last_sign:            # keep the larger of two same-sign extrema (noise)
            if abs(ypk) > ak[-1]:
                tk[-1], ak[-1] = t[i] + s * (t[i + 1] - t[i]), abs(ypk)
            continue
        tk.append(t[i] + s * (t[i + 1] - t[i]))
        ak.append(abs(ypk))
        last_sign = sign
    return np.asarray(tk), np.asarray(ak), theta0


def fit_envelope(tk: np.ndarray, ak: np.ndarray, I: float, mB: float, a_min: float = 0.0):
    """Half-period map A_{k+1} = e^{-b T/(4I)} A_k - 2 tau_c/(m B) (small amplitudes).

    Linear in (slope, intercept). Returns b, tau_c and the mean period T.
    """
    keep = ak[:-1] > a_min
    x, y = ak[:-1][keep], ak[1:][keep]
    if x.size < 3:
        raise ValueError("too few extrema for the envelope fit")
    A = np.vstack([x, np.ones_like(x)]).T
    (slope, icpt), *_ = np.linalg.lstsq(A, y, rcond=None)
    T = 2.0 * float(np.mean(np.diff(tk)))
    b = -4.0 * I * math.log(slope) / T if 0 < slope < 1 else float("nan")
    tau_c = max(0.0, -icpt) * mB / 2.0
    return b, tau_c, T


def integrate_decay(t: np.ndarray, I: float, mB: float, b: float, tau_c: float,
                    theta0: float, th_init: float, w_init: float, w_eps: float | None = None):
    """Equation of motion with a smoothed Coulomb term tau_c tanh(w / w_eps), sampled at t."""
    from scipy.integrate import solve_ivp
    w0 = math.sqrt(mB / I)
    w_eps = w_eps or 1e-3 * w0

    def rhs(_, y):
        th, w = y
        return [w, (-mB * math.sin(th - theta0) - b * w - tau_c * math.tanh(w / w_eps)) / I]

    sol = solve_ivp(rhs, (t[0], t[-1]), [th_init, w_init], t_eval=t, method="LSODA",
                    rtol=1e-8, atol=1e-10, max_step=0.05 / (w0 / (2 * math.pi)))
    return sol.y[0]


def fit_decay(t: np.ndarray, th: np.ndarray, I: float, mB: float, b0: float, tau_c0: float):
    """Least-squares fit of the whole record. Free: b, tau_c, theta0, theta(0), omega(0)."""
    from scipy.optimize import least_squares
    th0_guess = float(np.median(th[-max(3, len(th) // 10):]))
    k = max(1, min(5, len(th) // 50))
    w_init = float((th[k] - th[0]) / (t[k] - t[0])) if k > 1 else 0.0
    scale_b = max(b0, 1e-12)
    scale_c = max(tau_c0, 1e-3 * mB * 1e-3)

    def resid(p):
        b, tc, th0, thi, wi = p[0] * scale_b, p[1] * scale_c, p[2], p[3], p[4]
        return integrate_decay(t, I, mB, b, tc, th0, thi, wi) - th

    p0 = [max(b0, 0.0) / scale_b, max(tau_c0, 0.0) / scale_c, th0_guess, float(th[0]), w_init]
    # Finite-difference step well above the integrator tolerance, or the
    # Jacobian is noise and the optimizer does not move.
    res = least_squares(resid, p0, bounds=([0, 0, -np.inf, -np.inf, -np.inf], np.inf),
                        diff_step=1e-4, x_scale="jac")
    b, tc = res.x[0] * scale_b, res.x[1] * scale_c
    J = res.jac
    dof = max(1, len(th) - len(p0))
    s2 = 2 * res.cost / dof
    try:
        cov = np.linalg.inv(J.T @ J) * s2
        err = np.sqrt(np.diag(cov))
        eb, ec = err[0] * scale_b, err[1] * scale_c
    except np.linalg.LinAlgError:
        eb = ec = float("nan")
    return dict(b=b, b_err=eb, tau_c=tc, tau_c_err=ec, theta0=float(res.x[2]),
                rms_rad=float(math.sqrt(s2)), success=bool(res.success))


def stuck_index(t: np.ndarray, th_smooth: np.ndarray, w0: float, frac: float = 0.01,
                noise: float = 0.0) -> int:
    """First sample after which the (smoothed) angle moves less than frac of the
    initial amplitude over two periods: the needle has stuck. The ODE fit uses
    the record up to there plus half a period."""
    amp = np.max(np.abs(th_smooth - th_smooth[-1]))
    T = 2 * math.pi / w0
    dt = t[1] - t[0]
    n_win = max(3, int(2 * T / dt))
    for i in range(len(th_smooth) - n_win):
        if np.ptp(th_smooth[i:i + n_win]) < max(frac * amp, 4.0 * noise):
            return min(len(th_smooth), i + int(0.5 * T / dt))
    return len(th_smooth)


def simulate(I: float, m: float, B: float, b: float, tau_c: float, tau_s: float,
             theta_init: float, t_end: float, fps: float, noise_rad: float = 0.0,
             theta0: float = 0.0, seed: int = 0):
    """Stick-slip free decay sampled at fps (fine semi-implicit Euler internally)."""
    rng = np.random.default_rng(seed)
    mB = m * B
    w0 = math.sqrt(mB / I)
    dt = (2 * math.pi / w0) / 2000.0
    n_out = int(t_end * fps) + 1
    t_out = np.arange(n_out) / fps
    out = np.empty(n_out)
    th, w, t, j = theta_init, 0.0, 0.0, 0
    stuck = False
    while j < n_out:
        while j < n_out and t_out[j] <= t + 1e-15:
            out[j] = th
            j += 1
        tau_m = -mB * math.sin(th - theta0)
        if stuck:
            if abs(tau_m) > tau_s:
                stuck = False
            else:
                t += dt
                continue
        acc = (tau_m - b * w - tau_c * math.copysign(1.0, w if w != 0 else tau_m)) / I
        w_new = w + acc * dt
        if w != 0 and w_new * w < 0 and abs(mB * math.sin(th - theta0)) <= tau_s:
            w_new, stuck = 0.0, True
        w = w_new
        th += w * dt
        t += dt
    return t_out, out + noise_rad * rng.standard_normal(n_out)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _read(path: str, cols: list[str]) -> list[np.ndarray]:
    import pandas as pd
    df = pd.read_csv(path)
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(f"{path}: missing columns {missing}; has {list(df.columns)}")
    return [df[c].to_numpy(float) for c in cols]


def _print(d: dict) -> None:
    print(json.dumps(d, indent=2))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("geometry", help="I from weighed mass and dimensions")
    g.add_argument("--mass", type=float, required=True, help="blade mass [kg] (without pivot)")
    g.add_argument("--length", type=float, required=True, help="tip-to-tip length L [m]")
    g.add_argument("--width", type=float, required=True, help="maximum width W [m]")
    g.add_argument("--pivot_mass", type=float, default=0.0)
    g.add_argument("--pivot_radius", type=float, default=0.0)
    g.add_argument("--density", type=float, default=7850.0, help="for the effective thickness [kg/m^3]")

    a = sub.add_parser("addmass", help="I from the frequency shift of a known added inertia")
    a.add_argument("--f0", type=float, required=True, help="frequency without the added inertia [Hz]")
    a.add_argument("--f1", type=float, required=True, help="frequency with it, same field [Hz]")
    a.add_argument("--dI", type=float, required=True, help="added inertia [kg m^2]")

    f = sub.add_parser("freq", help="m/I from omega^2 versus field (CSV columns B_T, f_Hz)")
    f.add_argument("csv")
    f.add_argument("--I", type=float, default=None, help="if given, also report m")

    d = sub.add_parser("dipole", help="m from the on-axis field (CSV columns r_m, B_T)")
    d.add_argument("csv")

    c = sub.add_parser("decay", help="b and tau_c from a free decay (CSV columns t, theta)")
    c.add_argument("csv")
    c.add_argument("--I", type=float, required=True)
    c.add_argument("--m", type=float, required=True)
    c.add_argument("--B", type=float, required=True, help="total horizontal field at the needle [T]")
    c.add_argument("--deg", action="store_true", help="theta column in degrees")
    c.add_argument("--t_col", default="t")
    c.add_argument("--theta_col", default="theta")
    c.add_argument("--no_ode", action="store_true", help="envelope fit only")

    s = sub.add_parser("deadband", help="tau_s from final rest angles (CSV column theta_final)")
    s.add_argument("csv")
    s.add_argument("--m", type=float, required=True)
    s.add_argument("--B", type=float, required=True)
    s.add_argument("--deg", action="store_true")

    m = sub.add_parser("simulate", help="write a synthetic free decay to CSV (t, theta)")
    m.add_argument("out")
    m.add_argument("--I", type=float, default=1.63e-10)
    m.add_argument("--m", type=float, default=4.9e-3)
    m.add_argument("--B", type=float, default=1e-4)
    m.add_argument("--b", type=float, default=2e-10)
    m.add_argument("--tau_c", type=float, default=5e-9)
    m.add_argument("--tau_s", type=float, default=8e-9)
    m.add_argument("--theta_init_deg", type=float, default=30.0)
    m.add_argument("--t_end", type=float, default=8.0)
    m.add_argument("--fps", type=float, default=240.0)
    m.add_argument("--noise_deg", type=float, default=0.2)
    m.add_argument("--seed", type=int, default=0)

    args = ap.parse_args(argv)

    if args.cmd == "geometry":
        I = inertia_rhombus(args.mass, args.length, args.width, args.pivot_mass, args.pivot_radius)
        thick = args.mass / (args.density * 0.5 * args.length * args.width)
        _print(dict(I=I, effective_thickness_m=thick))
    elif args.cmd == "addmass":
        _print(dict(I=inertia_added_mass(args.f0, args.f1, args.dI)))
    elif args.cmd == "freq":
        B, fz = _read(args.csv, ["B_T", "f_Hz"])
        k, B_res, rel = fit_freq_vs_field(B, fz)
        out = dict(m_over_I=k, B_residual_T=B_res, rel_rms=rel)
        if args.I:
            out["m"] = k * args.I
        _print(out)
    elif args.cmd == "dipole":
        r, B = _read(args.csv, ["r_m", "B_T"])
        mm, bg = fit_dipole(r, B)
        _print(dict(m=mm, B_background_T=bg))
    elif args.cmd == "decay":
        t, th = _read(args.csv, [args.t_col, args.theta_col])
        if args.deg:
            th = np.deg2rad(th)
        mB = args.m * args.B
        w0 = math.sqrt(mB / args.I)
        ths = smooth(t, th, 2 * math.pi / w0)
        noise = float(np.std(th - ths))
        n = stuck_index(t, ths, w0, noise=noise)
        t, th, ths = t[:n], th[:n], ths[:n]
        tk, ak, th0 = find_extrema(t, ths, min_sep=0.25 * 2 * math.pi / w0)
        b_env, tc_env, T = fit_envelope(tk, ak, args.I, mB, a_min=5.0 * noise)
        out = dict(f0_Hz=w0 / (2 * math.pi), period_measured_s=T, n_extrema=int(len(ak)),
                   noise_rad=noise, samples_used=int(n),
                   envelope=dict(b=b_env, tau_c=tc_env))
        if not args.no_ode:
            out["ode_fit"] = fit_decay(t, th, args.I, mB, b_env, tc_env)
            bb = out["ode_fit"]["b"]
        else:
            bb = b_env
        out["Q_single_needle"] = math.sqrt(mB * args.I) / bb if bb > 0 else float("inf")
        out["tau_c_over_mB"] = (out.get("ode_fit", out["envelope"])["tau_c"]) / mB
        _print(out)
    elif args.cmd == "deadband":
        (thf,) = _read(args.csv, ["theta_final"])
        if args.deg:
            thf = np.deg2rad(thf)
        ts, th0, half = tau_static_from_deadband(thf, args.m, args.B)
        _print(dict(tau_s=ts, theta0_rad=th0, half_band_rad=half, tau_s_over_mB=ts / (args.m * args.B)))
    elif args.cmd == "simulate":
        t, th = simulate(args.I, args.m, args.B, args.b, args.tau_c, args.tau_s,
                         math.radians(args.theta_init_deg), args.t_end, args.fps,
                         math.radians(args.noise_deg), seed=args.seed)
        import pandas as pd
        pd.DataFrame(dict(t=t, theta=th)).to_csv(args.out, index=False)
        print(f"wrote {args.out} ({len(t)} samples)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
