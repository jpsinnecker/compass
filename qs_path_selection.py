#!/usr/bin/env python3
"""
qs_path_selection.py -- inertial path selection in the quasi-static loops.

In the event-driven quasi-static grid (fase2_campaign.py --protocol qs) the
loop of the gradient flow ("od", Q0 -> 0) and of the inertial dynamics at
Q0 = 10 coincide in most conditions, but differ by 10-43 % at a few
geometry/size combinations. This script finds where the two loops part.

Both runs of a pair start from the same random state (same seed) and are
identical at rest at +B_max (end of the first branch "0_to_pos"; saturation
erases the initial state, seeds 0 and 1 give the same loop). From there the
rows (one per field step, at rest) are aligned one to one, since the field
path is the same. With d_i = K (M_proj^Q10 - M_proj^od) at row i, in
"needle equivalents" (one full reversal changes K M_proj by 2):

  onset       first row with |d_i| > --onset_thr (states no longer equal;
              can be a smooth drift of canted needles, without reversals);
  divergence  first row with |d_i| > --div_thr (default 1 = half a reversal);
  jumps       rows where |d_i - d_{i-1}| > --jump_thr (a discrete event: the
              avalanche of that field step has a different outcome).

For each, the field, branch, avalanche size (flip_field, flip_angle) and
dissipated energy (E_diss increment / m B_ref) of both runs are reported.

Subcommands
  scan   align the existing od / Q10 CSVs, write path_selection_scan.csv and
         path_selection_jumps.csv;
  rerun  rerun selected pairs with --qs_state_log and --event_log (exact
         campaign arguments, so the CSVs are reproduced), optionally stopping
         after --stop_after field steps;
  map    from the reruns, map the rest states of od and Q10 on the lattice at
         the onset, the divergence and the largest jump; classify the needles
         that differ (edge vs bulk) and count multiple reversals within one
         field step (overshoot) from the event log. Writes
         path_selection_spatial.csv and figures.

    python3 qs_path_selection.py scan  --out_dir ~/compass_runs/fase2_qs_paths
    python3 qs_path_selection.py rerun --out_dir ~/compass_runs/fase2_qs_paths \\
        --cases honeycomb_N10 honeycomb_N06 --stop_after 0
    python3 qs_path_selection.py map   --out_dir ~/compass_runs/fase2_qs_paths
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent
DEFAULT_RUN_DIRS = [Path("~/compass_runs/fase2_qs").expanduser(),
                    Path("~/compass_runs/fase2_qs_mumax").expanduser()]
DEFAULT_OUT = Path("~/compass_runs/fase2_qs_paths").expanduser()
PAIR = ("od", "Q10")


# ---------------------------------------------------------------- loading

def find_run(tag: str, run_dirs) -> tuple[Path, Path] | None:
    for d in run_dirs:
        c, m = d / "data" / f"{tag}.csv", d / "meta" / f"{tag}.json"
        if c.exists() and m.exists():
            return c, m
    return None


def case_tag(case: str, label: str, seed: int) -> str:
    return f"{case}_qs_{label}_s{seed}"


def parse_case(case: str) -> tuple[str, int]:
    g, n = case.rsplit("_N", 1)
    return g, int(n)


def all_cases(run_dirs) -> list[str]:
    cases = set()
    for d in run_dirs:
        for p in (d / "data").glob("*_qs_od_s*.csv"):
            if not p.name.endswith("_events.csv"):
                cases.add(p.name.split("_qs_")[0])
    return sorted(cases)


def load_pair(case: str, seed: int, run_dirs):
    out = {}
    for lab in PAIR:
        r = find_run(case_tag(case, lab, seed), run_dirs)
        if r is None:
            return None
        df = pd.read_csv(r[0])
        meta = json.loads(r[1].read_text())
        out[lab] = (df, meta)
    return out


# ---------------------------------------------------------------- scan

def step_quantities(df: pd.DataFrame, meta: dict) -> pd.DataFrame:
    """Per-row avalanche observables: s (flip_field), s_angle, e = dE_diss / (m B_ref)."""
    norm = meta["config"]["moment"] * meta["derived"]["B_ref_T"]
    e = np.diff(df.E_diss_J.values, prepend=df.E_diss_J.values[0]) / norm
    return pd.DataFrame(dict(s=df.flip_field.values, s_angle=df.flip_angle.values, e=e))


def scan_pair(case: str, seed: int, run_dirs, onset_thr: float, div_thr: float, jump_thr: float):
    pair = load_pair(case, seed, run_dirs)
    if pair is None:
        return None, None
    (da, ma), (db, mb) = pair["od"], pair["Q10"]
    n = min(len(da), len(db))
    da, db = da.iloc[:n].reset_index(drop=True), db.iloc[:n].reset_index(drop=True)
    if not np.allclose(da.B_scalar_T.values, db.B_scalar_T.values, rtol=0, atol=1e-12):
        raise SystemExit(f"{case}: field paths differ between od and Q10")
    K = ma["derived"]["K"]
    B_ref = ma["derived"]["B_ref_T"]
    qa, qb = step_quantities(da, ma), step_quantities(db, mb)
    n_up = int((da.branch == "0_to_pos").sum())          # row n_up = rest at +B_max
    d = K * (db.M_proj.values - da.M_proj.values)
    jump = np.diff(d, prepend=d[0])
    rows = np.arange(n)
    after = rows > n_up

    def first(mask):
        idx = np.where(mask & after)[0]
        return int(idx[0]) if len(idx) else None

    def describe(i, kind):
        return dict(case=case, seed=seed, kind=kind, row=i, branch=da.branch[i],
                    B_over_Bref=da.B_scalar_T[i] / B_ref, K=K, d=d[i], jump=jump[i],
                    M_od=da.M_proj[i], M_Q10=db.M_proj[i],
                    s_od=qa.s[i], s_Q10=qb.s[i], s_angle_od=qa.s_angle[i], s_angle_Q10=qb.s_angle[i],
                    e_od=qa.e[i], e_Q10=qb.e[i])

    summary = dict(case=case, seed=seed, K=K, n_rows=n, row_Bmax=n_up, d_at_Bmax=d[n_up])
    jumps = []
    i_on, i_div = first(np.abs(d) > onset_thr), first(np.abs(d) > div_thr)
    i_jump = first(np.abs(jump) > jump_thr)
    for key, i in (("onset", i_on), ("divergence", i_div), ("first_jump", i_jump)):
        summary[f"row_{key}"] = i
        if i is not None:
            for k, v in describe(i, key).items():
                if k not in ("case", "seed", "kind", "row", "K"):
                    summary[f"{key}_{k}"] = v
    for i in np.where((np.abs(jump) > jump_thr) & after)[0]:
        jumps.append(describe(int(i), "jump"))
    if jumps:
        big = max(jumps, key=lambda r: abs(r["jump"]))
        summary["row_largest_jump"] = big["row"]
        summary["largest_jump"] = big["jump"]
        summary["largest_jump_B_over_Bref"] = big["B_over_Bref"]
    summary["n_jumps"] = len(jumps)
    summary["d_final"] = d[-1]
    # Seed / machine sensitivity: the two seeds are the same loop after +B_max
    # up to the rest tolerance (and some were run on another machine), so
    # their difference measures how sensitive each path is to tiny changes.
    other = 1 - seed
    for lab in PAIR:
        a = find_run(case_tag(case, lab, seed), run_dirs)
        b = find_run(case_tag(case, lab, other), run_dirs)
        if a is None or b is None:
            continue
        ma_, mb_ = pd.read_csv(a[0]).M_proj.values, pd.read_csv(b[0]).M_proj.values
        m = min(len(ma_), len(mb_))
        dd = K * np.abs(ma_[n_up:m] - mb_[n_up:m])
        summary[f"seed_maxdiff_{lab}"] = float(dd.max())
        big = np.where(dd > div_thr)[0]
        summary[f"seed_first_row_{lab}"] = int(big[0] + n_up) if len(big) else None
    return summary, jumps


def cmd_scan(args):
    run_dirs = [Path(p).expanduser() for p in args.run_dirs]
    out = Path(args.out_dir).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    cases = args.cases or all_cases(run_dirs)
    summ, jumps = [], []
    for c in cases:
        s, j = scan_pair(c, args.seed, run_dirs, args.onset_thr, args.div_thr, args.jump_thr)
        if s is None:
            continue
        summ.append(s)
        jumps += j
    summ = pd.DataFrame(summ)
    jumps = pd.DataFrame(jumps)
    summ.to_csv(out / "path_selection_scan.csv", index=False)
    jumps.to_csv(out / "path_selection_jumps.csv", index=False)
    cols = ["case", "row_onset", "row_divergence", "divergence_B_over_Bref", "divergence_branch",
            "row_first_jump", "first_jump_B_over_Bref", "first_jump_jump", "first_jump_s_od",
            "first_jump_s_Q10", "first_jump_e_od", "first_jump_e_Q10", "n_jumps", "largest_jump", "d_final",
            "seed_maxdiff_od", "seed_maxdiff_Q10"]
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.float_format", "{:.3g}".format):
        print(summ[[c for c in cols if c in summ.columns]].to_string(index=False))
    print(f"\nwrote {out / 'path_selection_scan.csv'} and {out / 'path_selection_jumps.csv'}")


# ---------------------------------------------------------------- rerun

def campaign_job(tag: str, run_dirs) -> dict:
    """The campaign job record (exact damping, dynamics) from the run manifests."""
    for d in run_dirs:
        man = d / "manifest.jsonl"
        if man.exists():
            for line in open(man):
                job = json.loads(line)
                if job["tag"] == tag and job.get("status") == "ok":
                    return job
    raise SystemExit(f"{tag}: not found in the manifests of {run_dirs}")


def rerun_cmd(case: str, label: str, seed: int, out_dir: Path, stop_after: int, run_dirs) -> list[str]:
    """Same arguments as fase2_campaign.run_job (so the CSV is reproduced), plus state/event logs."""
    import fase2_campaign as fc
    job = campaign_job(case_tag(case, label, seed), run_dirs)
    cmd = [
        sys.executable, str(REPO / "compass.py"),
        "--geometry", job["geometry"], "--N", str(job["N"]), "--M", str(job["N"]),
        "--B_max_factor", str(fc.B_MAX_FACTOR),
        "--needle_thickness", str(fc.NEEDLE_THICKNESS),
        "--damping", repr(job["damping"]), "--damping_noise", "0",
        "--field_mode", "hysteresis_qs", "--qs_dynamics", job["dynamics"],
        "--seed", str(seed), "--log_every", "1",
        "--out_dir", str(out_dir), "--tag", job["tag"], "--png_dpi", "30",
        "--qs_state_log", "--event_log",
    ]
    if stop_after > 0:
        cmd += ["--qs_stop_after", str(stop_after)]
    return cmd


def cmd_rerun(args):
    out = Path(args.out_dir).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    run_dirs = [Path(p).expanduser() for p in args.run_dirs]
    env = dict(os.environ, OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    jobs = [(c, lab) for c in args.cases for lab in PAIR]

    def run(job):
        c, lab = job
        tag = case_tag(c, lab, args.seed)
        if (out / "meta" / f"{tag}.json").exists() and not args.force:
            return f"{tag}: exists"
        cmd = rerun_cmd(c, lab, args.seed, out, args.stop_after, run_dirs)
        with open(out / f"{tag}.log", "w") as log:
            rc = subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, cwd=str(REPO)).returncode
        return f"{tag}: rc={rc}"

    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        for msg in ex.map(run, jobs):
            print(msg, flush=True)


# ---------------------------------------------------------------- map

def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def coordination(xs, ys, r_nn, tol=1.05):
    dx = xs[:, None] - xs[None, :]
    dy = ys[:, None] - ys[None, :]
    r = np.hypot(dx, dy)
    return ((r > 1e-12) & (r < tol * r_nn)).sum(axis=1)


def load_rerun(out: Path, case: str, label: str, seed: int):
    tag = case_tag(case, label, seed)
    st = np.load(out / "states" / f"{tag}_qs_states.npz")
    df = pd.read_csv(out / "data" / f"{tag}.csv")
    ev_p = out / "data" / f"{tag}_events.csv"
    ev = pd.read_csv(ev_p) if ev_p.exists() else pd.DataFrame(columns=["step", "t_s", "needle_id", "channel", "theta_rad"])
    meta = json.loads((out / "meta" / f"{tag}.json").read_text())
    return dict(xs=st["xs"], ys=st["ys"], r_nn=float(st["r_nn"]), theta=st["theta"],
                phi=float(st["phi_ext"]), df=df, ev=ev, meta=meta)


def check_reproduction(case, seed, run, run_dirs, label):
    """Max |M_proj| difference between a rerun and the campaign CSV (same rows)."""
    r = find_run(case_tag(case, label, seed), run_dirs)
    if r is None:
        return float("nan")
    ref = pd.read_csv(r[0])
    n = min(len(ref), len(run["df"]))
    return float(np.max(np.abs(ref.M_proj.values[:n] - run["df"].M_proj.values[:n])))


def step_events(run, row: int) -> pd.DataFrame:
    """Committed events of field step `row` (integration steps in (step_{row-1}, step_row])."""
    st = run["df"].step.values
    lo, hi = st[row - 1], st[row]
    ev = run["ev"]
    return ev[(ev.step > lo) & (ev.step <= hi)]


# Needle states are compared through the projection on the drive axis, as in
# the flip_field channel: a needle "reversed" (between two rest states, or
# between od and Q10) if sign(cos(theta - phi)) changed and it rotated by more
# than MIN_ROT_DEG (so that jitter of needles near the perpendicular does not
# count). Canted needles in kagome-like arrays reverse with rotations < 90 deg.
MIN_ROT_DEG = 30.0
CANT_DEG = 10.0       # "canted differently" above this


def reversed_between(th_a, th_b, phi):
    rot = np.abs(wrap(th_b - th_a)) > np.radians(MIN_ROT_DEG)
    return rot & (np.sign(np.cos(th_a - phi)) != np.sign(np.cos(th_b - phi)))


def timeline(ro, rq, coord) -> pd.DataFrame:
    """Per-row comparison of the od and Q10 rest states (rows both reruns have)."""
    n = min(len(ro["theta"]), len(rq["theta"]), len(ro["df"]), len(rq["df"]))
    th_o, th_q = ro["theta"][:n], rq["theta"][:n]
    K = th_o.shape[1]
    edge = coord < coord.max()
    phi = ro["phi"]
    dth = np.abs(wrap(th_q - th_o))
    mo = np.zeros_like(dth, dtype=bool)
    mq = np.zeros_like(dth, dtype=bool)
    mo[1:] = reversed_between(th_o[:-1], th_o[1:], phi)
    mq[1:] = reversed_between(th_q[:-1], th_q[1:], phi)
    differ = reversed_between(th_o, th_q, phi)
    B_ref = ro["meta"]["derived"]["B_ref_T"]
    norm = ro["meta"]["config"]["moment"] * B_ref
    tl = pd.DataFrame(dict(
        row=np.arange(n), branch=ro["df"].branch.values[:n],
        B_over_Bref=ro["df"].B_scalar_T.values[:n] / B_ref,
        d=K * (rq["df"].M_proj.values[:n] - ro["df"].M_proj.values[:n]),
        n_differ=differ.sum(1), n_differ_edge=(differ & edge).sum(1),
        n_cant=(dth > np.radians(CANT_DEG)).sum(1), max_dtheta_deg=np.degrees(dth.max(1)),
        moved_od=mo.sum(1), moved_od_edge=(mo & edge).sum(1),
        moved_Q10=mq.sum(1), moved_Q10_edge=(mq & edge).sum(1),
        s_od=ro["df"].flip_field.values[:n], s_Q10=rq["df"].flip_field.values[:n],
        e_od=np.diff(ro["df"].E_diss_J.values[:n], prepend=0.0) / norm,
        e_Q10=np.diff(rq["df"].E_diss_J.values[:n], prepend=0.0) / norm,
    ))
    # needles reversed more than once (drive-axis channel) within one field step
    for lab, run in (("od", ro), ("Q10", rq)):
        ev = run["ev"]
        ev = ev[ev.channel == "field"]
        rows = np.searchsorted(run["df"].step.values[:n], ev.step.values, side="left")
        cnt = pd.DataFrame(dict(row=rows, needle=ev.needle_id.values)).groupby(["row", "needle"]).size()
        multi = (cnt >= 2).groupby(level=0).sum()
        tl[f"multi_flip_{lab}"] = tl.row.map(multi).fillna(0).astype(int)
        tl[f"max_flips_{lab}"] = tl.row.map(cnt.groupby(level=0).max()).fillna(0).astype(int)
    return tl


# ---------------------------------------------------------------- stability

class Landscape:
    """Energy landscape of one array: E = -1/2 sum m_i.B_i^dip - sum m_i.B_ext."""

    def __init__(self, run):
        import compass as cs
        be = cs.select_backend(False)
        self.T = cs.precompute_dipolar_tensor(run["xs"], run["ys"], np.inf, 0.0, be)
        meta = run["meta"]
        self.m = meta["config"]["moment"]
        self.b = meta["config"]["damping"]
        self.B_ref = meta["derived"]["B_ref_T"]
        self.B_loc = meta["derived"]["B_loc_bound_T"]
        self.tol_tau = meta["quasi_static"]["tol_torque_Nm"]
        self.tol_om = meta["quasi_static"]["tol_omega_rad_s"]

    def field(self, th, bx, by):
        c, s = np.cos(th), np.sin(th)
        T, m = self.T, self.m
        return m * (T.Axx @ c + T.Axy @ s) + bx, m * (T.Axy @ c + T.Ayy @ s) + by

    def torque(self, th, bx, by):
        Bx, By = self.field(th, bx, by)
        return self.m * (np.cos(th) * By - np.sin(th) * Bx)

    def hessian(self, th, bx, by):
        """d2E/dtheta_i dtheta_j in units of m B_ref."""
        c, s = np.cos(th), np.sin(th)
        T, m = self.T, self.m
        Bx, By = self.field(th, bx, by)
        H = -m * m * (np.outer(s, s) * T.Axx - np.outer(s, c) * T.Axy - np.outer(c, s) * T.Axy.T
                      + np.outer(c, c) * T.Ayy)
        np.fill_diagonal(H, m * (c * Bx + s * By))
        return H / (m * self.B_ref)

    def low_modes(self, th, bx, by, k=4):
        from scipy.linalg import eigh
        H = self.hessian(th, bx, by)
        k = min(k, len(th))
        w = eigh(H, eigvals_only=True, subset_by_index=[0, k - 1])
        return w

    def relax_od(self, th, bx, by, cfl=0.2, check=5, max_steps=2_000_000):
        """Gradient flow b dtheta/dt = tau with the engine's explicit Euler step and rest criterion."""
        dt = cfl * self.b / (self.m * self.B_loc)
        th = th.copy()
        for k in range(1, max_steps + 1):
            tau = self.torque(th, bx, by)
            th = wrap(th + dt * tau / self.b)
            if k % check == 0:
                ta = np.abs(self.torque(th, bx, by)).max()
                if ta <= self.tol_tau and ta / self.b <= self.tol_om:
                    return th, k
        return th, -1


def stability_columns(tl, ro, rq, rows, k=4):
    """Lowest Hessian eigenvalue and number of unstable modes of both rest states."""
    land = Landscape(ro)
    out = {f"{c}_{lab}": np.full(len(tl), np.nan) for c in ("lam_min", "n_neg") for lab in PAIR}
    for r in rows:
        for lab, run in (("od", ro), ("Q10", rq)):
            df = run["df"]
            w = land.low_modes(run["theta"][r], df.Bx_T[r], df.By_T[r], k)
            out[f"lam_min_{lab}"][r] = w[0]
            # count of negative modes, saturating at k (only k lowest computed)
            out[f"n_neg_{lab}"][r] = int((w < -1e-9).sum())
    for c, v in out.items():
        tl[c] = v
    return tl, land


def noisy_od_control(land, ro, rq, r0, r1, eps, seed=0):
    """Gradient flow over field steps r0..r1 from the od rest state at r0-1,
    with a random kick of eps rad on every needle before each relaxation.
    A rest state that is a saddle is then left, as in any real array.
    Returns per-row: needles reversed relative to Q10 and to od, energies."""
    rng = np.random.default_rng(seed)
    th = ro["theta"][r0 - 1].copy()
    rows = []
    phi = ro["phi"]
    for r in range(r0, r1 + 1):
        bx, by = ro["df"].Bx_T[r], ro["df"].By_T[r]
        th, nsteps = land.relax_od(th + eps * rng.standard_normal(th.size), bx, by)
        w = land.low_modes(th, bx, by, 2)
        rows.append(dict(row=r, B_over_Bref=ro["df"].B_scalar_T[r] / land.B_ref,
                         differ_vs_Q10=int(reversed_between(th, rq["theta"][r], phi).sum()),
                         differ_vs_od=int(reversed_between(th, ro["theta"][r], phi).sum()),
                         d_vs_od=float(np.cos(th - phi).sum() - np.cos(ro["theta"][r] - phi).sum()),
                         d_Q10_vs_od=float(np.cos(rq["theta"][r] - phi).sum() - np.cos(ro["theta"][r] - phi).sum()),
                         lam_min=float(w[0]), relax_steps=nsteps))
    return pd.DataFrame(rows)


def key_events(tl: pd.DataFrame, row_Bmax: int, onset: int | None) -> list[tuple[str, int]]:
    """Rows worth mapping: first split (first avalanche that leaves needles in a
    different state), maximum split, reconvergence (if the states meet again)."""
    t = tl[tl.row > row_Bmax]
    ev = []
    if onset is not None and np.isfinite(onset):
        ev.append(("onset", int(onset)))
    split = t[(t.n_differ > 0) & ((t.moved_od > 0) | (t.moved_Q10 > 0))]
    if len(split):
        r0 = int(split.row.iloc[0])
        ev.append(("first_split", r0))
        after = t[(t.row > r0) & (t.n_differ == 0) & (t.n_cant == 0)]
        r1 = int(after.row.iloc[0]) if len(after) else int(t.row.iloc[-1])
        win = t[(t.row >= r0) & (t.row <= r1)]
        ev.append(("max_split", int(win.row[win.n_differ.idxmax()])))
        if len(after):
            ev.append(("reconvergence", r1))
    return ev


def plot_pair(fig_path, case, ro, rq, row, kind, rec):
    import matplotlib.pyplot as plt
    phi = ro["phi"]
    xs, ys = ro["xs"] / ro["r_nn"], ro["ys"] / ro["r_nn"]
    dth = np.degrees(np.abs(wrap(rq["theta"][row] - ro["theta"][row])))
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 4.3))
    L = 0.42
    for ax, run, title in ((axs[0], ro, r"gradient flow ($Q_0\to0$)"), (axs[1], rq, r"$Q_0=10$")):
        th, prev = run["theta"][row], run["theta"][row - 1]
        moved = reversed_between(prev, th, phi)
        ax.scatter(xs[moved], ys[moved], s=110, c="#f2c14e", zorder=0, linewidths=0,
                   label=f"reversed in this step ({moved.sum()})")
        ax.quiver(xs - 0.5 * L * np.cos(th), ys - 0.5 * L * np.sin(th), L * np.cos(th), L * np.sin(th),
                  np.cos(th - phi), cmap="coolwarm_r", clim=(-1, 1), angles="xy", scale_units="xy",
                  scale=1, width=0.005, headwidth=3.5, headlength=4, zorder=2)
        ax.set_title(title, fontsize=10)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.01), fontsize=8, frameon=False)
    sc = axs[2].scatter(xs, ys, c=dth, cmap="magma_r", vmin=0, vmax=180, s=28, edgecolors="#999999", linewidths=0.3)
    axs[2].set_title(r"$|\theta_{Q10}-\theta_{od}|$ [deg]", fontsize=10)
    fig.colorbar(sc, ax=axs[2], fraction=0.046, pad=0.02)
    for ax in axs:
        ax.set_aspect("equal")
        ax.set_xticks([]), ax.set_yticks([])
    g, n = parse_case(case)
    fig.suptitle(f"{g} {n}x{n} (K={rec['K']}), {kind}: field step {row}, branch {rec['branch']}, "
                 f"B/B_ref = {rec['B_over_Bref']:+.3f}; needles reversed in one run only: {rec['n_differ']} "
                 f"(edge {rec['n_differ_edge']}); field along +x (red = along, blue = against)", fontsize=9)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(fig_path.with_suffix(f".{ext}"), dpi=150)
    plt.close(fig)


def plot_timeline(fig_path, case, tl, row_Bmax, events):
    import matplotlib.pyplot as plt
    t = tl[tl.row > row_Bmax]
    fig, axs = plt.subplots(4, 1, figsize=(8, 9), sharex=True)
    x = t.row.values
    axs[0].plot(x, t.d, color="#0d366b", lw=1)
    axs[0].set_ylabel(r"$K\,(M^{Q10}-M^{od})$")
    axs[1].plot(x, t.n_differ, color="#c0392b", lw=1, label="reversed in one run only")
    axs[1].plot(x, t.n_cant, color="#e59866", lw=1, label=r"canted differently ($>10^\circ$)")
    axs[1].set_ylabel("needles")
    axs[1].legend(fontsize=8, frameon=False)
    axs[2].vlines(x, 0, t.moved_od, color="#5598e7", lw=1.2, label="od")
    axs[2].vlines(x, 0, -t.moved_Q10, color="#0d366b", lw=1.2, label="Q10 (down)")
    axs[2].axhline(0, color="k", lw=0.5)
    axs[2].set_ylabel("reversed in step")
    axs[2].legend(fontsize=8, frameon=False)
    if "lam_min_od" in t:
        for lab, col in (("od", "#9ec5f4"), ("Q10", "#0d366b")):
            ok = np.isfinite(t[f"lam_min_{lab}"].values)
            axs[3].plot(x[ok], t[f"lam_min_{lab}"].values[ok], ".-", ms=2, lw=0.8, color=col, label=lab)
        axs[3].axhline(0, color="k", lw=0.5)
        axs[3].set_yscale("symlog", linthresh=1e-2)
        axs[3].set_ylabel(r"$\lambda_{min}$ of Hessian [$mB_{ref}$]")
        axs[3].legend(fontsize=8, frameon=False)
    axs[-1].set_xlabel("field step (row)")
    for j, (kind, r) in enumerate(events):
        for ax in axs:
            ax.axvline(r, color="#888888", lw=0.6, ls=":")
        b = float(tl.B_over_Bref[tl.row == r].iloc[0])
        axs[0].annotate(f"{kind} ({b:+.2f})", (r, 0.97 - 0.12 * j), xycoords=("data", "axes fraction"),
                        fontsize=7, va="top", ha="left")
    fig.suptitle(f"{case}: od vs Q10 rest states per field step (labels: event, B/B_ref)", fontsize=10)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(fig_path.with_suffix(f".{ext}"), dpi=150)
    plt.close(fig)


def cmd_map(args):
    import matplotlib
    matplotlib.use("Agg")
    out = Path(args.out_dir).expanduser()
    run_dirs = [Path(p).expanduser() for p in args.run_dirs]
    figs = out / "figs"
    figs.mkdir(parents=True, exist_ok=True)
    scan = pd.read_csv(out / "path_selection_scan.csv")
    recs = []
    for case in args.cases:
        try:
            ro, rq = (load_rerun(out, case, lab, args.seed) for lab in PAIR)
        except FileNotFoundError:
            print(f"{case}: no rerun with --qs_state_log, skipped")
            continue
        rep = {lab: check_reproduction(case, args.seed, r, run_dirs, lab) for lab, r in zip(PAIR, (ro, rq))}
        print(f"{case}: reproduction of the campaign CSV, max|dM_proj| od={rep['od']:.2e} Q10={rep['Q10']:.2e}")
        sc = scan[(scan.case == case) & (scan.seed == args.seed)].iloc[0]
        row_Bmax = int(sc.row_Bmax)
        coord = coordination(ro["xs"], ro["ys"], ro["r_nn"])
        tl = timeline(ro, rq, coord)
        # onset from the rerun itself (it may differ from the campaign CSV if not reproduced)
        t = tl[tl.row > row_Bmax]
        on = t[np.abs(t.d) > args.onset_thr]
        events = key_events(tl, row_Bmax, int(on.row.iloc[0]) if len(on) else None)
        ev_rows = dict(events)
        K = ro["theta"].shape[1]
        stride = args.stability_stride or (1 if K <= 300 else 5)
        srows = sorted(set(range(row_Bmax, len(tl), stride)) | set(ev_rows.values())
                       | set(r - 1 for r in ev_rows.values()))
        tl, land = stability_columns(tl, ro, rq, srows)
        tl.to_csv(out / f"timeline_{case}.csv", index=False)
        if "onset" in ev_rows and args.control_eps:
            r0 = ev_rows["onset"] - 2
            r1 = ev_rows.get("reconvergence", ev_rows.get("max_split", r0) + 10)
            r1 = min(r1, len(tl) - 1)
            for eps in args.control_eps:
                ctl = noisy_od_control(land, ro, rq, r0, r1, eps)
                ctl.insert(0, "eps_rad", eps)
                ctl.insert(0, "case", case)
                ctl.to_csv(out / f"control_noisy_od_{case}_eps{eps:g}.csv", index=False)
                fin = ctl.iloc[-1]
                mid = ctl.loc[ctl.differ_vs_Q10.add(0).idxmax()]
                print(f"  noisy od control eps={eps:g} rad, rows {r0}-{r1}: max reversed vs od "
                      f"{ctl.differ_vs_od.max()}, vs Q10 {ctl.differ_vs_Q10.max()} (od vs Q10 max "
                      f"{int(tl.n_differ[r0:r1 + 1].max())}); unstable rest states: "
                      f"{int((ctl.lam_min < -1e-9).sum())}/{len(ctl)}")
        plot_timeline(figs / f"timeline_{case}", case, tl, row_Bmax, events)
        for kind, row in events:
            r = tl.loc[row].to_dict()
            rec = dict(case=case, seed=args.seed, kind=kind, K=ro["theta"].shape[1],
                       n_edge_sites=int((coord < coord.max()).sum()), repro_od=rep["od"], repro_Q10=rep["Q10"], **r)
            rec["relax_steps_od"] = int(ro["df"].step[row] - ro["df"].step[row - 1])
            rec["relax_steps_Q10"] = int(rq["df"].step[row] - rq["df"].step[row - 1])
            recs.append(rec)
            if kind != "onset" or args.map_onset:
                plot_pair(figs / f"map_{case}_{kind}_row{row}", case, ro, rq, row, kind, rec)
    df = pd.DataFrame(recs)
    sp = out / "path_selection_spatial.csv"
    if sp.exists() and len(df):            # keep the cases not redone in this call
        old = pd.read_csv(sp)
        df = pd.concat([old[~old.case.isin(df.case.unique())], df], ignore_index=True)
    df.to_csv(sp, index=False)
    cols = ["case", "kind", "row", "branch", "B_over_Bref", "K", "d", "n_differ", "n_differ_edge", "n_cant",
            "moved_od", "moved_od_edge", "moved_Q10", "moved_Q10_edge", "s_od", "s_Q10", "e_od", "e_Q10",
            "multi_flip_od", "multi_flip_Q10", "max_flips_Q10", "lam_min_od", "n_neg_od", "lam_min_Q10",
            "n_neg_Q10", "repro_Q10"]
    with pd.option_context("display.width", 250, "display.max_columns", 40, "display.float_format", "{:.3g}".format):
        print(df[cols].to_string(index=False))
    print(f"\nwrote {out / 'path_selection_spatial.csv'}, timeline_<case>.csv and figures in {figs}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--run_dirs", nargs="+", default=[str(p) for p in DEFAULT_RUN_DIRS])
    common.add_argument("--out_dir", default=str(DEFAULT_OUT))
    common.add_argument("--seed", type=int, default=0)
    common.add_argument("--cases", nargs="+", default=None, help="e.g. honeycomb_N10 triangular_N32 (default: all)")
    s = sub.add_parser("scan", parents=[common])
    s.add_argument("--onset_thr", type=float, default=0.1, help="|K dM| threshold for the onset [needle equivalents]")
    s.add_argument("--div_thr", type=float, default=1.0, help="|K dM| threshold for the divergence (1 = half a reversal)")
    s.add_argument("--jump_thr", type=float, default=0.5, help="|change of K dM| in one field step for a discrete event")
    s.set_defaults(fn=cmd_scan)
    r = sub.add_parser("rerun", parents=[common])
    r.add_argument("--stop_after", type=int, default=0, help="field steps to run (0 = full loop)")
    r.add_argument("--jobs", type=int, default=4)
    r.add_argument("--force", action="store_true")
    r.set_defaults(fn=cmd_rerun)
    m = sub.add_parser("map", parents=[common])
    m.add_argument("--onset_thr", type=float, default=0.1)
    m.add_argument("--stability_stride", type=int, default=0, help="rows between Hessian checks (0: 1 for K<=300, else 5)")
    m.add_argument("--control_eps", type=float, nargs="*", default=[1e-3], help="kick amplitudes [rad] of the noisy gradient-flow control (none to skip)")
    m.add_argument("--map_onset", action="store_true", help="also map the onset row (usually no reversal yet)")
    m.set_defaults(fn=cmd_map)
    args = ap.parse_args()
    if args.cmd in ("rerun", "map") and not args.cases:
        ap.error(f"{args.cmd} needs --cases")
    args.fn(args)


if __name__ == "__main__":
    main()
