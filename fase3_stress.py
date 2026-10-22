#!/usr/bin/env python3
"""
fase3_stress.py -- return-point-memory stress test of Phase 3
(design: docs/FASE3_PILOTO_RPM.md, section 7).

The pilot (fase3_pilot.py) found exact RPM for one turning point per
geometry. This test tries to break it: per geometry, 7-8 turning points A_i
spread over the whole descending reversal (M from about +0.75 to -0.75 in the
Phase-2 QS2 loops, including the fields where inertia selected different
paths), each followed by an excursion A_i -> A_i + D -> A_i, and a nested
minor loop at the middle point (A -> A+D -> A+D/2 -> A+D -> A):

    4 -> A_1 -> A_1+D -> A_1 -> A_2 -> ... -> A_k -> A_k+D -> A_k -> -4

Every return to the turning point where a minor loop started is compared
with the state there before the loop (rpm.rpm_returns(match="last")).

Grid: N = 6, 10, 16, 24 (gradient flow, Q0 = 3, 10) and N = 32 (gradient flow, Q0 = 10),
seeds 0 and 1 (the pilot showed no seed dependence of the RPM outcome).

  python3 fase3_stress.py run --out_dir ~/compass_runs/fase3_stress --workers 8 [--dry_run]
  python3 fase3_stress.py analyze --run_dirs ~/compass_runs/fase3_stress
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO))
import fase2_campaign as f2  # noqa: E402
import fase3_pilot as p3  # noqa: E402

# Turning points A_i (decreasing), excursion D and index of the nested loop,
# in units of B_ref, from the descending branches of the 24x24 and 32x32 QS2
# loops (fields where M crosses +0.75 ... -0.75); multiples of 0.005.
STRESS = {
    "square": ([0.36, 0.22, 0.10, -0.04, -0.18, -0.30, -0.45], 0.30, 3),
    "triangular": ([0.025, 0.01, -0.01, -0.04, -0.10, -0.20, -0.32, -0.37], 0.20, 3),
    "honeycomb": ([0.57, 0.04, -0.105, -0.16, -0.215, -0.25, -0.30, -0.58], 0.15, 3),
}
_ALL = [("od", 0.0), ("Q3", 3.0), ("Q10", 10.0)]
GRID = [(6, _ALL), (10, _ALL), (16, _ALL), (24, _ALL), (32, [("od", 0.0), ("Q10", 10.0)])]
SEEDS = [0, 1]
STEP_FRACTION = 0.85        # field steps relative to a full five-branch loop (cost estimate only)


def program_points(geometry: str) -> list[float]:
    A, D, nested = STRESS[geometry]
    pts = [p3.B_SAT]
    for i, a in enumerate(A):
        pts += [a, a + D, a + D / 2, a + D, a] if i == nested else [a, a + D, a]
    pts.append(-p3.B_SAT)
    pts = [round(p, 6) for p in pts]
    if any(x == y for x, y in zip(pts, pts[1:])):
        raise ValueError(f"{geometry}: repeated consecutive turning point")
    return pts


def program_text(geometry: str) -> str:
    return ",".join(f"{p:g}" for p in program_points(geometry))


def build_jobs(out_dir: Path, seeds=None) -> list[dict]:
    seeds = seeds if seeds is not None else SEEDS
    jobs = []
    for n, dyn in GRID:
        for g in STRESS:
            p = f2.lattice_params(g, n)
            for (label, q0), seed in ((d, s) for d in dyn for s in seeds):
                q_for_b = q0 if q0 > 0 else f2.QS_OD_Q0_FOR_TIME_UNIT
                steps = STEP_FRACTION * f2.QS_STEPS_N6[label] * (n / 6.0) ** f2.QS_SIZE_EXP[label]
                jobs.append(dict(
                    tag=f"{g}_N{n:02d}_{label}_s{seed}_stress", geometry=g, N=n, K=p["K"], label=label,
                    Q0_target=q0, dynamics="overdamped" if label == "od" else "inertial", seed=seed,
                    kind="stress", program=program_text(g),
                    damping=math.sqrt(p["moment"] * p["B_lat"] * p["inertia"]) / q_for_b,
                    cost_est_s=steps * f2.us_per_step(n, False) * 1e-6, out_dir=str(out_dir)))
    jobs.sort(key=lambda j: (-j["cost_est_s"], j["tag"]))
    return jobs


def cmd_run(args) -> int:
    out_dir = Path(args.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    jobs = build_jobs(out_dir, seeds=[int(x) for x in args.seeds.split(",")] if args.seeds else None)
    if args.max_N:
        jobs = [j for j in jobs if j["N"] <= args.max_N]
    total = sum(j["cost_est_s"] for j in jobs)
    print(f"{len(jobs)} jobs; rough CPU estimate {total/3600:.1f} h (M2 Max per-step cost; "
          f"Phase-2 QS runs took ~2x this under full load)")
    for g in STRESS:
        print(f"  {g:10s} {program_text(g)}")
    (out_dir / "stress_plan.json").write_text(json.dumps(dict(
        programs={g: program_text(g) for g in STRESS}, grid=GRID, seeds=SEEDS, n_jobs=len(jobs),
        cli=vars(args), git_commit=subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                                                  capture_output=True, text=True).stdout.strip()), indent=2))
    if args.dry_run:
        return 0
    manifest = out_dir / "manifest.jsonl"
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(p3.run_job, j) for j in jobs]
        for done, fut in enumerate(as_completed(futures), start=1):
            res = fut.result()
            with open(manifest, "a") as fh:
                fh.write(json.dumps(res) + "\n")
            print(f"[{done}/{len(jobs)}] {time.time()-t0:7.0f}s  {res['status']:>16s}  {res['tag']}", flush=True)
    return 0


def cmd_analyze(args) -> int:
    import numpy as np
    import pandas as pd
    import rpm

    rows, per_return = [], []
    for rd in args.run_dirs:
        rd = os.path.expanduser(rd)
        for meta in sorted(Path(rd, "meta").glob("*_stress.json")):
            run = rpm.load_run(rd, meta.stem)
            g, n, label, seed = meta.stem.split("_")[:4]
            ret = rpm.rpm_returns(run, match="last")
            tp = rpm.turning_points(run)
            # Needles reversed during each minor loop (start state vs its far end):
            # a return after a loop that reversed nothing is a trivial test.
            far = {}
            for k in range(1, len(tp) - 1):
                far[int(tp.leg[k])] = rpm.state_distance(run.theta[tp.row[k]], run.theta[tp.row[k - 1]],
                                                         run.phi_ext)["flip_mismatch"]
            for _, r in ret.iterrows():
                per_return.append(dict(geometry=g, N=int(n[1:]), label=label, seed=int(seed[1:]),
                                       B_Bref=r.B_Bref, return_leg=int(r.return_leg),
                                       flipped_in_loop=far.get(int(r.return_leg) - 1, np.nan),
                                       flip_mismatch=r.flip_mismatch, rms_rad=r.rms_rad, dM=r.dM))
            rows.append(dict(geometry=g, N=int(n[1:]), label=label, seed=int(seed[1:]), returns=len(ret),
                             violating=int((ret.flip_mismatch > 0).sum()),
                             max_flip_mismatch=float(ret.flip_mismatch.max()),
                             max_rms_rad=float(ret.rms_rad.max()),
                             unconverged=run.meta["quasi_static"]["unconverged_field_steps"]))
    if not rows:
        print("no stress runs found")
        return 1
    out = Path(os.path.expanduser(args.out or args.run_dirs[0]))
    df, pr = pd.DataFrame(rows), pd.DataFrame(per_return)
    df.to_csv(out / "rpm_stress_runs.csv", index=False)
    pr.to_csv(out / "rpm_stress_returns.csv", index=False)
    pd.set_option("display.width", 200)
    print(df.sort_values(["N", "geometry", "label", "seed"]).round(5).to_string(index=False))
    nontriv = pr[pr.flipped_in_loop > 0]
    print(f"\n{len(pr)} returns tested, {len(nontriv)} after a minor loop that reversed needles; "
          f"{int((pr.flip_mismatch > 0).sum())} with needles in a different state")
    print(f"written {out/'rpm_stress_runs.csv'} and {out/'rpm_stress_returns.csv'}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--out_dir", required=True)
    r.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    r.add_argument("--dry_run", action="store_true")
    r.add_argument("--seeds", default=None)
    r.add_argument("--max_N", type=int, default=None, help="restrict to N <= max_N")
    a = sub.add_parser("analyze")
    a.add_argument("--run_dirs", nargs="+", required=True)
    a.add_argument("--out", default=None)
    args = ap.parse_args(argv)
    return cmd_run(args) if args.cmd == "run" else cmd_analyze(args)


if __name__ == "__main__":
    raise SystemExit(main())
