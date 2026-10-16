#!/usr/bin/env python3
"""
fase2_campaign.py -- Phase 2 of docs/AUDIT_PHYSICS_ROADMAP.md: dynamic
hysteresis map A(r, Q0) for the three geometries and several array sizes.

For every (geometry, N, Q0, r, seed) it runs one full five-branch hysteresis
cycle of compass.py at the dimensionless sweep rate r (--sweep_rate_dimless).
The damping b is chosen so that the engine's *intrinsic* Q0 (metadata "Q",
based on the lattice field B_lat of that geometry and size) equals the target
exactly: b = sqrt(m B_lat I) / Q0.

Each job is an independent `python3 compass.py ...` subprocess with
single-threaded BLAS, so --workers processes use that many cores. Jobs are
ordered longest-first for load balancing. A job is skipped if its metadata
file already exists (resume after an interruption). One JSON line per finished
job is appended to <out_dir>/manifest.jsonl.

Analysis: fase2_analysis.py.

Example:
    python3 fase2_campaign.py --out_dir ~/compass_runs/fase2 --workers 10
    python3 fase2_campaign.py --out_dir ~/compass_runs/fase2 --dry_run
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

import numpy as np

REPO = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO))
import compass as cs  # noqa: E402

GEOMETRIES = ["square", "triangular", "honeycomb"]
SIZES = [6, 10, 16, 24, 32]
Q0_VALUES = [0.1, 0.3, 1.0, 3.0, 10.0]
R_VALUES = [1e-1, 3e-2, 1e-2, 3e-3, 1e-3]
R_EXTRA_SMALL = [3e-4]          # only for N <= SMALL_N_MAX
SMALL_N_MAX = 16
SEEDS = [0]
SEED_CHECK = {"N": 10, "seeds": [1, 2]}   # seed-robustness subset
B_MAX_FACTOR = 4.0
NEEDLE_THICKNESS = 0.0004       # m; passed explicitly (see USER_GUIDE, thickness note)
TARGET_ROWS = 4000              # approximate CSV rows per run

# Calibrated on an M2 Max (triangular, B_max = 4 B_ref, r = 0.1, Q0 ~ 0.6):
# ~4500 steps; steps scale as (0.1 / r) * max(1, 1/Q0) / 1.6.
STEPS_AT_R01 = 4500.0
US_PER_STEP = {6: 47, 10: 56, 16: 103, 24: 300, 32: 800}


def lattice_params(geometry: str, N: int) -> dict:
    """m, I, B_lat, B_ref for a geometry/size, computed exactly as compass.py does."""
    geom = cs.make_lattice(N, N, geometry, cs.R_DEFAULT, 0.8, 0.0, np.random.default_rng(0))
    moment = cs.compute_moment_from_geometry(geom.needle_len, geom.needle_width, NEEDLE_THICKNESS)
    inertia = cs.compute_inertia_from_geometry(geom.needle_len, geom.needle_width, NEEDLE_THICKNESS)
    backend = cs.select_backend(False)
    tensor = cs.precompute_dipolar_tensor(geom.xs, geom.ys, math.inf, math.inf, backend)
    B_lat = moment * cs.lattice_field_per_moment(tensor, backend)
    B_ref = cs.MU0_OVER_4PI * 2.0 * moment / geom.r_nn ** 3
    return dict(K=geom.K, moment=moment, inertia=inertia, B_lat=B_lat, B_ref=B_ref)


def build_jobs(out_dir: Path) -> list[dict]:
    params = {(g, n): lattice_params(g, n) for g in GEOMETRIES for n in SIZES}
    jobs = []
    for g, n, q0 in itertools.product(GEOMETRIES, SIZES, Q0_VALUES):
        p = params[(g, n)]
        damping = math.sqrt(p["moment"] * p["B_lat"] * p["inertia"]) / q0
        r_list = R_VALUES + (R_EXTRA_SMALL if n <= SMALL_N_MAX else [])
        seeds = SEEDS + (SEED_CHECK["seeds"] if n == SEED_CHECK["N"] else [])
        for r, seed in itertools.product(r_list, seeds):
            steps = STEPS_AT_R01 * (0.1 / r) * max(1.0, 1.0 / q0) / 1.6
            tag = f"{g}_N{n:02d}_Q{q0:g}_r{r:g}_s{seed}"
            jobs.append(dict(
                tag=tag, geometry=g, N=n, K=p["K"], Q0_target=q0, r=r, seed=seed,
                damping=damping, B_lat=p["B_lat"], B_ref=p["B_ref"],
                steps_est=steps, cost_est_s=steps * US_PER_STEP[n] * 1e-6,
                log_every=max(1, int(steps / TARGET_ROWS)),
                out_dir=str(out_dir),
            ))
    jobs.sort(key=lambda j: -j["cost_est_s"])
    return jobs


def run_job(job: dict) -> dict:
    meta = Path(job["out_dir"]) / "meta" / f"{job['tag']}.json"
    if meta.exists():
        return dict(job, status="skipped_existing")
    env = dict(os.environ, OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    cmd = [
        sys.executable, str(REPO / "compass.py"),
        "--geometry", job["geometry"], "--N", str(job["N"]), "--M", str(job["N"]),
        "--field_mode", "hysteresis", "--B_max_factor", str(B_MAX_FACTOR),
        "--needle_thickness", str(NEEDLE_THICKNESS),
        "--damping", repr(job["damping"]), "--damping_noise", "0",
        "--sweep_rate_dimless", repr(job["r"]),
        "--seed", str(job["seed"]), "--log_every", str(job["log_every"]),
        "--out_dir", job["out_dir"], "--tag", job["tag"], "--png_dpi", "30",
    ]
    t0 = time.time()
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True, cwd=str(REPO))
    status = "ok" if proc.returncode == 0 else f"error:{proc.returncode}"
    out = dict(job, status=status, wall_s=time.time() - t0)
    if proc.returncode != 0:
        out["stderr_tail"] = proc.stderr[-2000:]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--dry_run", action="store_true", help="print the plan and the cost estimate only")
    ap.add_argument("--max_N", type=int, default=None, help="restrict to N <= max_N (pilot runs)")
    ap.add_argument("--only_r", type=str, default=None, help="comma list restricting r values (pilot runs)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    jobs = build_jobs(out_dir)
    if args.max_N is not None:
        jobs = [j for j in jobs if j["N"] <= args.max_N]
    if args.only_r:
        keep = {float(x) for x in args.only_r.split(",")}
        jobs = [j for j in jobs if j["r"] in keep]

    total = sum(j["cost_est_s"] for j in jobs)
    longest = max(j["cost_est_s"] for j in jobs)
    print(f"{len(jobs)} jobs; estimated CPU time {total/3600:.1f} h; "
          f"wall >= max({total/args.workers/3600:.2f} h, longest job {longest/60:.0f} min) "
          f"with {args.workers} workers")
    (out_dir / "campaign_plan.json").write_text(json.dumps(dict(
        geometries=GEOMETRIES, sizes=SIZES, Q0=Q0_VALUES, r=R_VALUES, r_extra_small=R_EXTRA_SMALL,
        small_N_max=SMALL_N_MAX, seeds=SEEDS, seed_check=SEED_CHECK, B_max_factor=B_MAX_FACTOR,
        needle_thickness=NEEDLE_THICKNESS, n_jobs=len(jobs), est_cpu_h=total / 3600,
        git_commit=subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True,
                                  text=True).stdout.strip(),
    ), indent=2))
    if args.dry_run:
        return 0

    manifest = out_dir / "manifest.jsonl"
    done = 0
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run_job, j) for j in jobs]
        for fut in as_completed(futures):
            res = fut.result()
            done += 1
            with open(manifest, "a") as fh:
                fh.write(json.dumps(res) + "\n")
            print(f"[{done}/{len(jobs)}] {time.time()-t0:7.0f}s  {res['status']:>16s}  {res['tag']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
