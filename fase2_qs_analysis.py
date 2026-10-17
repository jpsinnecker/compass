#!/usr/bin/env python3
"""
fase2_qs_analysis.py -- analysis of the event-driven quasi-static grid
(fase2_campaign.py --protocol qs): quasi-static loop area a0 per geometry,
size, relaxation dynamics (gradient flow "od", inertial Q0 = 1, 3, 10) and
seed; inertial effect relative to the gradient flow; per-field-step
avalanche statistics.

Per run (one CSV row per field step at rest):
  a      loop area per needle [m B_ref] from the energy ledger over the closed
         cycle (from the rest state at +B_max, end of "0_to_pos", to the end);
  Bc, Mr coercive field / B_ref and remanence (fase2_analysis.run_metrics);
  events per field step: s = flip_field count of the row (needles that
         reversed with respect to the drive axis), e = E_diss increment of
         the row / (m B_ref) (energy dissipated by the step).

Outputs in <analysis_dir>: qs_runs.csv, qs_summary.csv, qs_events.csv and
PDF/PNG figures in <analysis_dir>/figs (and copies for the report with
--report_figs docs/figs).

    python3 fase2_qs_analysis.py --run_dirs ~/compass_runs/fase2_qs ~/compass_runs/fase2_qs_mumax \\
        --analysis_dir ~/compass_runs/fase2_qs_analysis --report_figs docs/figs
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

import fase2_analysis as fa

GEOM_ORDER = fa.GEOM_ORDER
LABELS = ["od", "Q1", "Q3", "Q10"]
LABEL_TEXT = {"od": r"$Q_0\to 0$ (gradient flow)", "Q1": r"$Q_0=1$", "Q3": r"$Q_0=3$", "Q10": r"$Q_0=10$"}
LABEL_COLOR = {"od": "#9ec5f4", "Q1": "#5598e7", "Q3": "#256abf", "Q10": "#0d366b"}


def collect(run_dirs):
    runs, events, seen = [], [], set()
    for d in run_dirs:
        man = d / "manifest.jsonl"
        if not man.exists():
            continue
        for line in open(man):
            job = json.loads(line)
            tag = job["tag"]
            if job.get("status") not in ("ok", "skipped_existing") or tag in seen:
                continue
            meta_p, csv_p = d / "meta" / f"{tag}.json", d / "data" / f"{tag}.csv"
            if not (meta_p.exists() and csv_p.exists()):
                continue
            seen.add(tag)
            meta = json.loads(meta_p.read_text())
            m = fa.run_metrics(csv_p, meta)
            label = job.get("qs_label") or tag.split("_qs_")[1].rsplit("_s", 1)[0]
            runs.append(dict(tag=tag, geometry=job["geometry"], N=job["N"], K=job["K"], label=label,
                             seed=job["seed"], run_dir=str(d), **m))
            df = pd.read_csv(csv_p)
            norm = meta["config"]["moment"] * meta["derived"]["B_ref_T"]
            e = np.diff(df.E_diss_J.values, prepend=df.E_diss_J.values[0]) / norm
            ev = pd.DataFrame(dict(tag=tag, geometry=job["geometry"], N=job["N"], label=label,
                                   seed=job["seed"], branch=df.branch.values,
                                   B_over_Bref=df.B_scalar_T.values / meta["derived"]["B_ref_T"],
                                   s=df.flip_field.values, s_angle=df.flip_angle.values, e=e))
            events.append(ev.iloc[1:])          # row 0 is the initial state
    runs = pd.DataFrame(runs).sort_values(["geometry", "N", "label", "seed"]).reset_index(drop=True)
    events = pd.concat(events, ignore_index=True) if events else pd.DataFrame()
    return runs, events


def summarize(runs):
    g = runs.groupby(["geometry", "N", "label"])
    s = g.agg(a_mean=("a", "mean"), a_min=("a", "min"), a_max=("a", "max"), n_seeds=("a", "count"),
              Bc_mean=("Bc", "mean"), Mr_mean=("Mr", "mean"),
              unconverged=("qs_unconverged", "sum"), max_residual=("balance_residual_rel", "max")).reset_index()
    s["seed_spread"] = (s.a_max - s.a_min) / s.a_mean
    od = s[s.label == "od"].set_index(["geometry", "N"]).a_mean
    s["rel_to_od"] = [r.a_mean / od.get((r.geometry, r.N), np.nan) - 1.0 for r in s.itertuples()]
    return s


def _style(ax):
    fa._style(ax)


def figures(summ, events, fig_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Fig Q1: a0 vs N per geometry, one line per dynamics (seed range as bars).
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, g in zip(axes, GEOM_ORDER):
        for lab in LABELS:
            d = summ[(summ.geometry == g) & (summ.label == lab)].sort_values("N")
            if d.empty:
                continue
            ax.errorbar(d.N, d.a_mean, yerr=[d.a_mean - d.a_min, d.a_max - d.a_mean], fmt="-o",
                        color=LABEL_COLOR[lab], lw=2, ms=4, capsize=2, label=LABEL_TEXT[lab])
        ax.set_title(g, fontsize=10, color="#222222")
        ax.set_xlabel("array size $N$", fontsize=9)
        _style(ax)
    axes[0].set_ylabel(r"quasi-static loop area $a_0$ [$m B_{\rm ref}$]", fontsize=9)
    axes[-1].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(fig_dir / f"figQS1_a0_vs_N.{ext}", dpi=200)
    plt.close(fig)

    # Fig Q2: inertial effect a0(Q0)/a0(od) - 1 vs N.
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
    for ax, g in zip(axes, GEOM_ORDER):
        ax.axhline(0, color="#999999", lw=0.8)
        for lab in LABELS[1:]:
            d = summ[(summ.geometry == g) & (summ.label == lab)].sort_values("N")
            if d.empty:
                continue
            ax.plot(d.N, 100 * d.rel_to_od, "-o", color=LABEL_COLOR[lab], lw=2, ms=4, label=LABEL_TEXT[lab])
        ax.set_title(g, fontsize=10, color="#222222")
        ax.set_xlabel("array size $N$", fontsize=9)
        _style(ax)
    axes[0].set_ylabel(r"$a_0(Q_0)/a_0(Q_0\to0)-1$  [%]", fontsize=9)
    axes[-1].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(fig_dir / f"figQS2_inertial_effect.{ext}", dpi=200)
    plt.close(fig)

    # Fig Q3: CCDF of per-step avalanche size s (> 0), largest N, per dynamics.
    if not events.empty:
        nmax = events.N.max()
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
        for ax, g in zip(axes, GEOM_ORDER):
            for lab in LABELS:
                s = events[(events.geometry == g) & (events.N == nmax) & (events.label == lab)].s.values
                s = np.sort(s[s > 0])
                if s.size == 0:
                    continue
                ccdf = 1.0 - np.arange(s.size) / s.size
                ax.step(s, ccdf, where="post", color=LABEL_COLOR[lab], lw=1.8, label=f"{LABEL_TEXT[lab]} (n={s.size})")
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_title(f"{g}, N={nmax}", fontsize=10, color="#222222")
            ax.set_xlabel("avalanche size $s$ (reversed needles per field step)", fontsize=8)
            _style(ax)
        axes[0].set_ylabel(r"$P(S \geq s)$", fontsize=9)
        for ax in axes:
            ax.legend(fontsize=6, frameon=False)
        fig.tight_layout()
        for ext in ("pdf", "png"):
            fig.savefig(fig_dir / f"figQS3_avalanche_ccdf.{ext}", dpi=200)
        plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run_dirs", nargs="+", required=True)
    ap.add_argument("--analysis_dir", required=True)
    ap.add_argument("--report_figs", default=None, help="also copy the PDF figures here (e.g. docs/figs)")
    args = ap.parse_args()
    an = Path(args.analysis_dir).expanduser()
    figs = an / "figs"
    figs.mkdir(parents=True, exist_ok=True)

    runs, events = collect([Path(d).expanduser() for d in args.run_dirs])
    summ = summarize(runs)
    runs.to_csv(an / "qs_runs.csv", index=False)
    summ.to_csv(an / "qs_summary.csv", index=False)
    events.to_csv(an / "qs_events.csv", index=False)
    print(f"{len(runs)} QS runs; unconverged field steps total = {int(runs.qs_unconverged.sum())}; "
          f"max balance residual = {runs.balance_residual_rel.max():.1e}")
    figures(summ, events, figs)
    if args.report_figs:
        dest = Path(args.report_figs)
        dest.mkdir(parents=True, exist_ok=True)
        for p in figs.glob("figQS*.pdf"):
            shutil.copy(p, dest / p.name)
    with pd.option_context("display.width", 170, "display.max_rows", 200):
        print(summ[["geometry", "N", "label", "n_seeds", "a_mean", "seed_spread", "rel_to_od",
                    "Bc_mean", "Mr_mean", "unconverged"]].round(4).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
