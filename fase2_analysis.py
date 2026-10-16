#!/usr/bin/env python3
"""
fase2_analysis.py -- analysis of the Phase 2 dynamic-hysteresis campaign
(fase2_campaign.py): loop metrics per run, quasi-static extrapolation
A0(Q0) = lim_{r->0} A(r, Q0), the rate exponent beta, and figures.

Per run (one closed cycle, from the first arrival at +B_max to the end, where
the field is +B_max again):
  a      loop area per needle in units of m*B_ref, from the energy ledger:
         a = [W_field(end) - W_field(t_top)] / (K m B_ref)  (W = field work);
  a_MdB  the same from the curve, -oint M_proj dB / B_ref (cross-check);
  d      dissipation per needle over the cycle, [E_diss(end)-E_diss(t_top)]/(K m B_ref);
  Bc     coercive field / B_ref (mean of |crossings| on the two branches);
  Mr     remanence (mean |M_proj| at B = 0 on the two descending-to-zero branches).

Per (geometry, N, Q0): fit a(r) = a0 + c * r**beta on all r (scipy curve_fit,
a0 >= 0, 0 < beta < 2), plus a model-free convergence flag: relative change of
a between the two smallest r.

Outputs in <out_dir>/analysis/: runs.csv, fits.csv, seed_check.csv and PNG/PDF
figures.

    python3 fase2_analysis.py --out_dir ~/compass_runs/fase2
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

GEOM_ORDER = ["square", "triangular", "honeycomb"]
# Reference palette (dataviz skill): categorical slots 1-3 for geometries,
# one-hue blue ramp (light -> dark) for the ordered Q0 values.
GEOM_COLOR = {"square": "#2a78d6", "triangular": "#eb6834", "honeycomb": "#1baf7a"}
GEOM_MARKER = {"square": "s", "triangular": "^", "honeycomb": "h"}
Q_RAMP = ["#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"]


# ----------------------------------------------------------------------------
# per-run metrics
# ----------------------------------------------------------------------------

def _crossing(B, M):
    """Field values where M changes sign (linear interpolation)."""
    s = np.sign(M)
    idx = np.where(s[:-1] * s[1:] < 0)[0]
    return [B[i] - M[i] * (B[i + 1] - B[i]) / (M[i + 1] - M[i]) for i in idx]


def run_metrics(csv_path: Path, meta: dict) -> dict:
    df = pd.read_csv(csv_path)
    der, sr = meta["derived"], meta["sweep_rate"]
    K, m, B_ref = der["K"], meta["config"]["moment"], der["B_ref_T"]
    Bmax = abs(meta["config"]["B_ext"])
    t_top = Bmax / sr["sweep_rate_T_per_s"]          # first arrival at +B_max
    t, W, D = df.t_s.values, df.W_field_J.values, df.E_diss_J.values
    W_top, D_top = np.interp(t_top, t, W), np.interp(t_top, t, D)
    norm = K * m * B_ref

    cyc = df[df.t_s >= t_top]
    B, M = cyc.B_scalar_T.values, cyc.M_proj.values
    a_MdB = -np.trapz(M, B) / B_ref

    desc = df[df.branch.isin(["pos_to_0", "0_to_neg"])]
    asc = df[df.branch.isin(["neg_to_0", "0_to_pos_final"])]
    c_desc, c_asc = _crossing(desc.B_scalar_T.values, desc.M_proj.values), _crossing(asc.B_scalar_T.values, asc.M_proj.values)
    Bc = np.nan
    if c_desc and c_asc:
        Bc = 0.5 * (abs(c_desc[-1]) + abs(c_asc[-1])) / B_ref
    Mr = []
    for br in ("pos_to_0", "neg_to_0"):
        seg = df[df.branch == br]
        if len(seg):
            Mr.append(abs(seg.M_proj.values[-1]))
    return dict(
        a=(W[-1] - W_top) / norm,
        a_MdB=a_MdB,
        d=(D[-1] - D_top) / norm,
        Bc=Bc,
        Mr=float(np.mean(Mr)) if Mr else np.nan,
        M_at_Bmax=float(df.M_proj.max()),
        Q0=der["Q"],
        r=sr["sweep_rate_dimless"],
        balance_residual_rel=meta["energy_ledger"]["balance_residual_rel"],
        flagged_steps=meta["stability_monitor"]["flagged_steps"],
        n_steps=der["n_steps"],
        runtime_s=meta.get("runtime_s", np.nan),
    )


def collect(out_dir: Path) -> pd.DataFrame:
    rows = []
    for line in open(out_dir / "manifest.jsonl"):
        job = json.loads(line)
        if job.get("status") not in ("ok", "skipped_existing"):
            continue
        meta_path = out_dir / "meta" / f"{job['tag']}.json"
        csv_path = out_dir / "data" / f"{job['tag']}.csv"
        if not (meta_path.exists() and csv_path.exists()):
            continue
        meta = json.loads(meta_path.read_text())
        rec = dict(tag=job["tag"], geometry=job["geometry"], N=job["N"], K=job["K"],
                   Q0_target=job["Q0_target"], r_target=job["r"], seed=job["seed"])
        rec.update(run_metrics(csv_path, meta))
        rows.append(rec)
    df = pd.DataFrame(rows).drop_duplicates("tag")
    return df.sort_values(["geometry", "N", "Q0_target", "r_target", "seed"]).reset_index(drop=True)


# ----------------------------------------------------------------------------
# quasi-static extrapolation
# ----------------------------------------------------------------------------

def _model(r, a0, c, beta):
    return a0 + c * r ** beta


def fit_condition(r, a):
    order = np.argsort(r)
    r, a = np.asarray(r)[order], np.asarray(a)[order]
    out = dict(n_r=len(r), r_min=r[0], a_at_rmin=a[0],
               rel_change_two_smallest=abs(a[1] - a[0]) / max(abs(a[0]), 1e-300) if len(r) > 1 else np.nan)
    if len(r) < 4:
        return dict(out, a0=np.nan, a0_err=np.nan, c=np.nan, beta=np.nan, beta_err=np.nan, fit_ok=False)
    try:
        p, cov = curve_fit(_model, r, a, p0=[max(a[0] * 0.8, 1e-6), (a[-1] - a[0]) / r[-1] ** 0.6, 0.6],
                           bounds=([0.0, -np.inf, 1e-3], [np.inf, np.inf, 2.0]), maxfev=20000)
        err = np.sqrt(np.diag(cov))
        return dict(out, a0=p[0], a0_err=err[0], c=p[1], beta=p[2], beta_err=err[2], fit_ok=True)
    except Exception:
        return dict(out, a0=np.nan, a0_err=np.nan, c=np.nan, beta=np.nan, beta_err=np.nan, fit_ok=False)


def fit_all(runs: pd.DataFrame) -> pd.DataFrame:
    base = runs[runs.seed == 0]
    rows = []
    for (g, n, q), grp in base.groupby(["geometry", "N", "Q0_target"]):
        f = fit_condition(grp.r.values, grp.a.values)
        fb = fit_condition(grp.r.values, grp.Bc.values)
        rows.append(dict(geometry=g, N=n, Q0=q, **f,
                         Bc0=fb["a0"], Bc_at_rmin=fb["a_at_rmin"], Bc_beta=fb["beta"]))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------------

def _style(ax):
    ax.grid(True, color="#e6e6e6", linewidth=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors="#444444", labelsize=8)


def figures(runs: pd.DataFrame, fits: pd.DataFrame, fig_dir: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    base = runs[runs.seed == 0]
    Qs = sorted(base.Q0_target.unique())
    qcol = {q: Q_RAMP[i % len(Q_RAMP)] for i, q in enumerate(Qs)}

    # Fig 1: a(r) per geometry at the largest N that has the smallest r.
    n_ref = max(n for n in base.N.unique() if base[base.N == n].r_target.min() <= 3.1e-4) \
        if (base.r_target <= 3.1e-4).any() else base.N.max()
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, g in zip(axes, GEOM_ORDER):
        for q in Qs:
            d = base[(base.geometry == g) & (base.N == n_ref) & (base.Q0_target == q)].sort_values("r")
            if d.empty:
                continue
            ax.plot(d.r, d.a, "-o", color=qcol[q], lw=2, ms=4, label=f"$Q_0$={q:g}")
            f = fits[(fits.geometry == g) & (fits.N == n_ref) & (fits.Q0 == q)]
            if len(f) and f.fit_ok.iloc[0]:
                ax.axhline(f.a0.iloc[0], color=qcol[q], lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_title(f"{g}, N={n_ref}", fontsize=10, color="#222222")
        ax.set_xlabel("sweep rate $r$", fontsize=9)
        _style(ax)
    axes[0].set_ylabel(r"loop area per needle $a$ [$m B_{\rm ref}$]", fontsize=9)
    axes[-1].legend(fontsize=7, frameon=False, title="dotted: $a_0$", title_fontsize=7)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig1_area_vs_rate.png", dpi=200)
    fig.savefig(fig_dir / "fig1_area_vs_rate.pdf")
    plt.close(fig)

    # Fig 2: quasi-static a0(Q0), one panel per geometry, all sizes (finite-size).
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    Ns = sorted(fits.N.unique())
    ncol = {n: Q_RAMP[min(i, len(Q_RAMP) - 1)] for i, n in enumerate(Ns)}
    for ax, g in zip(axes, GEOM_ORDER):
        for n in Ns:
            f = fits[(fits.geometry == g) & (fits.N == n)].sort_values("Q0")
            ok = f[f.fit_ok]
            ax.errorbar(ok.Q0, ok.a0, yerr=ok.a0_err, fmt="-o", color=ncol[n], lw=2, ms=4,
                        capsize=2, label=f"N={n}")
            ax.plot(f.Q0, f.a_at_rmin, "x", color=ncol[n], ms=5)
        ax.set_xscale("log")
        ax.set_title(g, fontsize=10, color="#222222")
        ax.set_xlabel("$Q_0$", fontsize=9)
        _style(ax)
    axes[0].set_ylabel(r"quasi-static loop area $a_0$ [$m B_{\rm ref}$]", fontsize=9)
    axes[-1].legend(fontsize=7, frameon=False, title="x: smallest r", title_fontsize=7)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig2_a0_vs_Q0.png", dpi=200)
    fig.savefig(fig_dir / "fig2_a0_vs_Q0.pdf")
    plt.close(fig)

    # Fig 3: excess area a - a0 vs r (log-log) -> rate exponent beta.
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, g in zip(axes, GEOM_ORDER):
        for q in Qs:
            d = base[(base.geometry == g) & (base.N == n_ref) & (base.Q0_target == q)].sort_values("r")
            f = fits[(fits.geometry == g) & (fits.N == n_ref) & (fits.Q0 == q)]
            if d.empty or not len(f) or not f.fit_ok.iloc[0]:
                continue
            ex = d.a.values - f.a0.iloc[0]
            m = ex > 0
            ax.plot(d.r.values[m], ex[m], "-o", color=qcol[q], lw=2, ms=4,
                    label=f"$Q_0$={q:g}, β={f.beta.iloc[0]:.2f}")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(f"{g}, N={n_ref}", fontsize=10, color="#222222")
        ax.set_xlabel("sweep rate $r$", fontsize=9)
        _style(ax)
    axes[0].set_ylabel(r"$a - a_0$", fontsize=9)
    for ax in axes:
        ax.legend(fontsize=6.5, frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig3_excess_area.png", dpi=200)
    fig.savefig(fig_dir / "fig3_excess_area.pdf")
    plt.close(fig)

    # Fig 4: example loops M(B) at Q0 = 0.1 and 10, three rates (square, n_ref).
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.6), sharey=True)
    return_figs = []
    for ax, q in zip(axes, [min(Qs), max(Qs)]):
        d = base[(base.geometry == "square") & (base.N == n_ref) & (base.Q0_target == q)].sort_values("r")
        rs = d.r_target.values
        pick = [rs[0], rs[len(rs) // 2], rs[-1]] if len(rs) >= 3 else rs
        for i, rr in enumerate(pick):
            tag = d[d.r_target == rr].tag.iloc[0]
            c = pd.read_csv(fig_dir.parent.parent / "data" / f"{tag}.csv")
            meta = json.loads((fig_dir.parent.parent / "meta" / f"{tag}.json").read_text())
            B_ref = meta["derived"]["B_ref_T"]
            sr = meta["sweep_rate"]
            t_top = abs(meta["config"]["B_ext"]) / sr["sweep_rate_T_per_s"]
            c = c[c.t_s >= t_top]
            ax.plot(c.B_scalar_T / B_ref, c.M_proj, color=Q_RAMP[[0, 2, 4][i]], lw=1.6, label=f"r={rr:g}")
        ax.set_title(f"square N={n_ref}, $Q_0$={q:g}", fontsize=10, color="#222222")
        ax.set_xlabel(r"$B/B_{\rm ref}$", fontsize=9)
        _style(ax)
        ax.legend(fontsize=7, frameon=False)
    axes[0].set_ylabel(r"$M_{\rm proj}$", fontsize=9)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig4_loops.png", dpi=200)
    fig.savefig(fig_dir / "fig4_loops.pdf")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()
    out_dir = Path(args.out_dir).expanduser()
    an = out_dir / "analysis"
    figs = an / "figs"
    figs.mkdir(parents=True, exist_ok=True)

    runs = collect(out_dir)
    runs.to_csv(an / "runs.csv", index=False)
    print(f"{len(runs)} runs collected")
    q = runs
    print(f"QC: max |Q0-target|/target = {np.nanmax(np.abs(q.Q0 / q.Q0_target - 1)):.1e}; "
          f"max balance residual = {q.balance_residual_rel.max():.1e}; flagged runs = {(q.flagged_steps > 0).sum()}; "
          f"min M at B_max = {q.M_at_Bmax.min():.3f}; max |a - a_MdB| = {np.nanmax(np.abs(q.a - q.a_MdB)):.2e}")

    fits = fit_all(runs)
    fits.to_csv(an / "fits.csv", index=False)

    seeds = runs[runs.N == runs[runs.seed > 0].N.min()] if (runs.seed > 0).any() else pd.DataFrame()
    if len(seeds):
        sc = seeds.groupby(["geometry", "Q0_target", "r_target"]).a.agg(["min", "max", "count"]).reset_index()
        sc["rel_spread"] = (sc["max"] - sc["min"]) / sc["max"].abs()
        sc.to_csv(an / "seed_check.csv", index=False)
        print(f"seed check: max relative spread of a across seeds = {sc.rel_spread.max():.2e}")

    figures(runs, fits, figs)
    with pd.option_context("display.width", 160, "display.max_rows", 200):
        print(fits[["geometry", "N", "Q0", "a_at_rmin", "a0", "a0_err", "beta", "beta_err",
                    "rel_change_two_smallest", "Bc_at_rmin"]].round(4).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
