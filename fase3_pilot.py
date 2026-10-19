#!/usr/bin/env python3
"""
fase3_pilot.py -- return-point-memory (RPM) pilot of Phase 3
(design: docs/FASE3_PILOTO_RPM.md).

Per geometry, a quasi-static field program with two minor excursions from a
turning point A in the middle of the descending reversal,

    4 -> A -> A+dp -> A -> A+dg -> A -> -4      (units of B_ref)

and a reference run 4 -> -4 with the same seed. Dynamics: gradient flow
("od", control without inertia) and Q0 = 10. All runs log the rest state of
every field step (--qs_state_log).

Tests (rpm.py), per geometry / dynamics / seed:
  RPM-1  state at A after the small excursion vs first passage (leg 4 vs 2)
  RPM-2  state at A after the large excursion vs first passage (leg 6 vs 2)
  RPM-3  outer branch after the excursions vs the reference run at equal
         fields below A (program leg 7 vs reference leg 2)

  python3 fase3_pilot.py run --out_dir ~/compass_runs/fase3_pilot --workers 7
  python3 fase3_pilot.py run --out_dir ... --dry_run
  python3 fase3_pilot.py analyze --run_dirs ~/compass_runs/fase3_pilot
"""
from __future__ import annotations

import argparse
import itertools
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

N_PILOT = 24
DB_FACTOR = 0.005
# A, small and large excursion (units of B_ref), from the 24x24 QS2 loops
# (docs/FASE3_PILOTO_RPM.md, section 2); all multiples of DB_FACTOR.
PROGRAMS = {
    "square": (-0.10, 0.145, 0.43),
    "triangular": (-0.20, 0.095, 0.28),
    "honeycomb": (-0.20, 0.075, 0.22),
}
B_SAT = 4.0
DYNAMICS = [("od", 0.0), ("Q10", 10.0)]
SEEDS = [0, 1, 2, 3, 4]
STEP_FRACTION = 0.6        # field steps relative to a full five-branch loop (cost estimate only)


def program_text(geometry: str) -> str:
    A, dp, dg = PROGRAMS[geometry]
    pts = [B_SAT, A, A + dp, A, A + dg, A, -B_SAT]
    return ",".join(f"{round(p, 6):g}" for p in pts)


def build_jobs(out_dir: Path, sizes=None, labels=None, seeds=None) -> list[dict]:
    sizes = sizes or [N_PILOT]
    seeds = seeds if seeds is not None else SEEDS
    dyn = [d for d in DYNAMICS if not labels or d[0] in labels]
    jobs = []
    for g, n in itertools.product(PROGRAMS, sizes):
        p = f2.lattice_params(g, n)
        for (label, q0), seed in itertools.product(dyn, seeds):
            q_for_b = q0 if q0 > 0 else f2.QS_OD_Q0_FOR_TIME_UNIT
            damping = math.sqrt(p["moment"] * p["B_lat"] * p["inertia"]) / q_for_b
            steps = STEP_FRACTION * f2.QS_STEPS_N6[label] * (n / 6.0) ** f2.QS_SIZE_EXP[label]
            for kind, prog in (("prog", program_text(g)), ("ref", f"{B_SAT:g},{-B_SAT:g}")):
                jobs.append(dict(
                    tag=f"{g}_N{n:02d}_{label}_s{seed}_{kind}", geometry=g, N=n, K=p["K"],
                    label=label, Q0_target=q0, dynamics="overdamped" if label == "od" else "inertial",
                    seed=seed, kind=kind, program=prog, damping=damping,
                    cost_est_s=steps * f2.us_per_step(n, False) * 1e-6, out_dir=str(out_dir)))
    jobs.sort(key=lambda j: (-j["cost_est_s"], j["tag"]))
    return jobs


def run_job(job: dict) -> dict:
    if (Path(job["out_dir"]) / "meta" / f"{job['tag']}.json").exists():
        return dict(job, status="skipped_existing")
    env = dict(os.environ, OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    cmd = [sys.executable, str(REPO / "compass.py"),
           "--geometry", job["geometry"], "--N", str(job["N"]), "--M", str(job["N"]),
           "--B_max_factor", str(B_SAT), "--needle_thickness", str(f2.NEEDLE_THICKNESS),
           "--damping", repr(job["damping"]), "--damping_noise", "0",
           "--field_mode", "hysteresis_qs", "--qs_dynamics", job["dynamics"],
           "--qs_dB_factor", str(DB_FACTOR), "--qs_program", job["program"], "--qs_state_log",
           "--seed", str(job["seed"]), "--log_every", "1",
           "--out_dir", job["out_dir"], "--tag", job["tag"], "--png_dpi", "30"]
    t0 = time.time()
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True, cwd=str(REPO))
    out = dict(job, status="ok" if proc.returncode == 0 else f"error:{proc.returncode}", wall_s=time.time() - t0)
    if proc.returncode != 0:
        out["stderr_tail"] = proc.stderr[-2000:]
    return out


def cmd_run(args) -> int:
    out_dir = Path(args.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    jobs = build_jobs(out_dir,
                      sizes=[int(x) for x in args.sizes.split(",")] if args.sizes else None,
                      labels=args.labels.split(",") if args.labels else None,
                      seeds=[int(x) for x in args.seeds.split(",")] if args.seeds else None)
    total = sum(j["cost_est_s"] for j in jobs)
    print(f"{len(jobs)} jobs; rough CPU estimate {total/3600:.1f} h (M2 Max per-step cost; "
          f"the Phase-2 QS runs took ~2x this under full load)")
    for g in PROGRAMS:
        print(f"  {g:10s} program {program_text(g)}")
    (out_dir / "pilot_plan.json").write_text(json.dumps(dict(
        programs={g: program_text(g) for g in PROGRAMS}, dynamics=DYNAMICS, seeds=SEEDS,
        dB_factor=DB_FACTOR, n_jobs=len(jobs), cli=vars(args),
        git_commit=subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True,
                                  text=True).stdout.strip()), indent=2))
    if args.dry_run:
        return 0
    manifest = out_dir / "manifest.jsonl"
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run_job, j) for j in jobs]
        for done, fut in enumerate(as_completed(futures), start=1):
            res = fut.result()
            with open(manifest, "a") as fh:
                fh.write(json.dumps(res) + "\n")
            print(f"[{done}/{len(jobs)}] {time.time()-t0:7.0f}s  {res['status']:>16s}  {res['tag']}", flush=True)
    return 0


def cmd_analyze(args) -> int:
    import pandas as pd
    import rpm

    found = {}
    for rd in args.run_dirs:
        rd = os.path.expanduser(rd)
        for meta in sorted(Path(rd, "meta").glob("*.json")):
            found[meta.stem] = rd
    rows = []
    for tag, rd in sorted(found.items()):
        if not tag.endswith("_prog") or tag[:-5] + "_ref" not in found:
            continue
        run = rpm.load_run(rd, tag)
        ref = rpm.load_run(found[tag[:-5] + "_ref"], tag[:-5] + "_ref")
        g, n, label, seed = tag.split("_")[:4]
        ret = rpm.rpm_returns(run)
        r1 = ret[ret.return_leg == 4].iloc[0]
        r2 = ret[ret.return_leg == 6].iloc[0]
        r3 = rpm.compare_legs(run, 7, ref, 2)
        rows.append(dict(
            geometry=g, N=int(n[1:]), label=label, seed=int(seed[1:]),
            rpm1_flip=r1.flip_mismatch, rpm1_rms=r1.rms_rad, rpm1_dM=r1.dM,
            rpm2_flip=r2.flip_mismatch, rpm2_rms=r2.rms_rad, rpm2_dM=r2.dM,
            rpm3_flip_max=float(r3.flip_mismatch.max()), rpm3_dM_maxabs=float(r3.dM.abs().max()),
            rpm3_n_fields=int(len(r3)),
            unconverged=run.meta["quasi_static"]["unconverged_field_steps"]
            + ref.meta["quasi_static"]["unconverged_field_steps"]))
    df = pd.DataFrame(rows)
    if df.empty:
        print("no complete program/reference pairs found")
        return 1
    out = Path(os.path.expanduser(args.out or args.run_dirs[0])) / "rpm_pilot.csv"
    df.to_csv(out, index=False)
    pd.set_option("display.width", 200)
    print(df.round(4).to_string(index=False))
    viol = df.assign(v1=df.rpm1_flip > 0, v2=df.rpm2_flip > 0, v3=df.rpm3_flip_max > 0)
    summ = viol.groupby(["geometry", "label"]).agg(
        seeds=("seed", "count"), rpm1_viol=("v1", "sum"), rpm2_viol=("v2", "sum"), rpm3_viol=("v3", "sum"),
        rpm2_flip_mean=("rpm2_flip", "mean"), rpm3_flip_mean=("rpm3_flip_max", "mean"))
    print("\nseeds violating each test (flip mismatch > 0):")
    print(summ.round(4).to_string())
    print(f"\nwritten {out}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--out_dir", required=True)
    r.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    r.add_argument("--dry_run", action="store_true")
    r.add_argument("--sizes", default=None, help="comma list of N (default 24)")
    r.add_argument("--labels", default=None, help="comma subset of od,Q10")
    r.add_argument("--seeds", default=None, help="comma list of seeds (default 0-4)")
    a = sub.add_parser("analyze")
    a.add_argument("--run_dirs", nargs="+", required=True)
    a.add_argument("--out", default=None, help="directory for rpm_pilot.csv (default: first run dir)")
    args = ap.parse_args(argv)
    return cmd_run(args) if args.cmd == "run" else cmd_analyze(args)


if __name__ == "__main__":
    raise SystemExit(main())
