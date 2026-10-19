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

Splitting a campaign between machines: every machine builds the same job
list and keeps only its shard (--shard k/n). Shards are balanced by estimated
cost with longest-processing-time assignment, weighted by --shard_weights
(relative throughput of each machine). Use the same --shard_weights and grid
options on every machine; run outputs can be merged afterwards by copying the
data/ and meta/ folders and concatenating the manifests.

Large arrays on a GPU (--use_gpu): runs one job at a time on the GPU
(--workers is forced to 1). Use it with --sizes for the finite-size extension,
e.g. --sizes 48,64,96.

Examples:
    python3 fase2_campaign.py --out_dir ~/compass_runs/fase2 --workers 10
    python3 fase2_campaign.py --out_dir ~/compass_runs/fase2 --dry_run
    # two machines, the first ~1.5x faster:
    python3 fase2_campaign.py --out_dir ~/compass_runs/fase2 --shard 1/2 --shard_weights 0.6,0.4 --workers 10
    python3 fase2_campaign.py --out_dir ~/compass_runs/fase2 --shard 2/2 --shard_weights 0.6,0.4 --workers 7
    # finite-size extension on the GPU:
    python3 fase2_campaign.py --out_dir ~/compass_runs/fase2_large --sizes 48,64,96 --use_gpu
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


def us_per_step(N: int, gpu: bool) -> float:
    """Rough per-step cost [us]. CPU: measured table (M2 Max), K^2 scaling beyond it.
    GPU: calibrated on an RTX 3090 (float64): 2900 us/step at K = 4096
    (square 64x64, 13 s run including CuPy kernel compilation), modelled as
    150 us of per-step overhead plus a K^2 term."""
    K = N * N
    if gpu:
        return 150.0 + 2750.0 * (K / 4096.0) ** 2
    if N in US_PER_STEP:
        return US_PER_STEP[N]
    return US_PER_STEP[32] * (K / 1024.0) ** 2


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


def build_jobs(out_dir: Path, sizes=None, q0_values=None, r_values=None, gpu: bool = False,
               seed_check: bool = True) -> list[dict]:
    sizes = sizes or SIZES
    q0_values = q0_values or Q0_VALUES
    r_base = r_values or R_VALUES
    params = {(g, n): lattice_params(g, n) for g in GEOMETRIES for n in sizes}
    jobs = []
    for g, n, q0 in itertools.product(GEOMETRIES, sizes, q0_values):
        p = params[(g, n)]
        damping = math.sqrt(p["moment"] * p["B_lat"] * p["inertia"]) / q0
        r_list = list(r_base) + ([] if r_values else (R_EXTRA_SMALL if n <= SMALL_N_MAX else []))
        seeds = SEEDS + (SEED_CHECK["seeds"] if (seed_check and n == SEED_CHECK["N"]) else [])
        for r, seed in itertools.product(r_list, seeds):
            steps = STEPS_AT_R01 * (0.1 / r) * max(1.0, 1.0 / q0) / 1.6
            tag = f"{g}_N{n:02d}_Q{q0:g}_r{r:g}_s{seed}"
            jobs.append(dict(
                tag=tag, geometry=g, N=n, K=p["K"], Q0_target=q0, r=r, seed=seed,
                damping=damping, B_lat=p["B_lat"], B_ref=p["B_ref"],
                steps_est=steps, cost_est_s=steps * us_per_step(n, gpu) * 1e-6, use_gpu=gpu,
                log_every=max(1, int(steps / TARGET_ROWS)),
                out_dir=str(out_dir),
            ))
    jobs.sort(key=lambda j: (-j["cost_est_s"], j["tag"]))
    return jobs


# Event-driven quasi-static protocol (--protocol qs): instead of the rate
# axis, the relaxation dynamics at fixed field: the gradient-flow limit
# (Q0 -> 0, "od") and inertial dynamics at several Q0. Steps per loop measured
# on 6x6 arrays (dB = 0.005 B_ref, tol = 1e-5); QS_SIZE_EXP scales them with N.
QS_DYNAMICS = [("od", 0.0), ("Q1", 1.0), ("Q3", 3.0), ("Q10", 10.0)]
# With the saddle-check kick (--qs_kick 1e-4), measured on honeycomb 6x6:
# od 2.8e6, Q3 3.6e6, Q10 11.4e6 steps; Q1 interpolated.
QS_STEPS_N6 = {"od": 2.8e6, "Q1": 2.6e6, "Q3": 3.6e6, "Q10": 11.4e6}
# Measured growth with N (square 6 -> 16): od ~N^0.46, Q3 ~N^0.1, Q10 ~N^0; Q1
# assumed in between until measured.
QS_SIZE_EXP = {"od": 0.46, "Q1": 0.3, "Q3": 0.1, "Q10": 0.0}
QS_SEEDS = [0, 1]
QS_OD_Q0_FOR_TIME_UNIT = 0.3     # damping used to set the time unit of the gradient flow


def build_qs_jobs(out_dir: Path, sizes=None, labels=None, seeds=None, gpu: bool = False) -> list[dict]:
    sizes = sizes or SIZES
    seeds = seeds or QS_SEEDS
    dyn = [d for d in QS_DYNAMICS if not labels or d[0] in labels]
    params = {(g, n): lattice_params(g, n) for g in GEOMETRIES for n in sizes}
    jobs = []
    for g, n, (label, q0), seed in itertools.product(GEOMETRIES, sizes, dyn, seeds):
        p = params[(g, n)]
        q_for_b = q0 if q0 > 0 else QS_OD_Q0_FOR_TIME_UNIT
        damping = math.sqrt(p["moment"] * p["B_lat"] * p["inertia"]) / q_for_b
        steps = QS_STEPS_N6[label] * (n / 6.0) ** QS_SIZE_EXP[label]
        jobs.append(dict(
            tag=f"{g}_N{n:02d}_qs_{label}_s{seed}", protocol="qs", geometry=g, N=n, K=p["K"],
            Q0_target=q0, dynamics="overdamped" if label == "od" else "inertial", qs_label=label,
            r=0.0, seed=seed, damping=damping, B_lat=p["B_lat"], B_ref=p["B_ref"],
            steps_est=steps, cost_est_s=steps * us_per_step(n, gpu) * 1e-6, use_gpu=gpu,
            log_every=1, out_dir=str(out_dir),
        ))
    jobs.sort(key=lambda j: (-j["cost_est_s"], j["tag"]))
    return jobs


def shard_jobs(jobs: list[dict], k: int, n: int, weights: list[float]) -> list[dict]:
    """Deterministic weighted LPT split: return the jobs of shard k (1-based)."""
    load = [0.0] * n
    mine = []
    for j in jobs:                       # already sorted by decreasing cost
        s = min(range(n), key=lambda i: ((load[i] + j["cost_est_s"]) / weights[i], i))
        load[s] += j["cost_est_s"]
        if s == k - 1:
            mine.append(j)
    return mine


def run_job(job: dict) -> dict:
    meta = Path(job["out_dir"]) / "meta" / f"{job['tag']}.json"
    if meta.exists():
        return dict(job, status="skipped_existing")
    env = dict(os.environ, OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    cmd = [
        sys.executable, str(REPO / "compass.py"),
        "--geometry", job["geometry"], "--N", str(job["N"]), "--M", str(job["N"]),
        "--B_max_factor", str(B_MAX_FACTOR),
        "--needle_thickness", str(NEEDLE_THICKNESS),
        "--damping", repr(job["damping"]), "--damping_noise", "0",
    ] + (["--field_mode", "hysteresis_qs", "--qs_dynamics", job["dynamics"]]
         if job.get("protocol") == "qs" else
         ["--field_mode", "hysteresis", "--sweep_rate_dimless", repr(job["r"])]) + [
        "--seed", str(job["seed"]), "--log_every", str(job["log_every"]),
        "--out_dir", job["out_dir"], "--tag", job["tag"], "--png_dpi", "30",
    ] + (["--use_gpu"] if job.get("use_gpu") else [])
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
    ap.add_argument("--sizes", type=str, default=None, help="comma list of N overriding the default sizes")
    ap.add_argument("--q0", type=str, default=None, help="comma list of Q0 overriding the default values")
    ap.add_argument("--rates", type=str, default=None, help="comma list of r overriding the default rates")
    ap.add_argument("--no_seed_check", action="store_true", help="omit the extra seeds of the seed-robustness subset")
    ap.add_argument("--use_gpu", action="store_true", help="run each job with compass.py --use_gpu, one at a time")
    ap.add_argument("--protocol", choices=["sweep", "qs"], default="sweep",
                    help="sweep: rate sweeps (r axis); qs: event-driven quasi-static loops (dynamics axis: od, Q1, Q3, Q10)")
    ap.add_argument("--qs_labels", type=str, default=None, help="comma subset of od,Q1,Q3,Q10 (qs protocol)")
    ap.add_argument("--qs_seeds", type=str, default=None, help="comma list of seeds (qs protocol; default 0,1)")
    ap.add_argument("--shard", type=str, default=None, help="k/n: run only shard k of n (1-based)")
    ap.add_argument("--shard_weights", type=str, default=None, help="comma list of n relative machine speeds (default equal)")
    ap.add_argument("--skip_tags", type=str, default=None, help="comma list of run tags already done on another machine (applied after sharding)")
    args = ap.parse_args()

    def floats(txt):
        return [float(x) for x in txt.split(",")] if txt else None

    out_dir = Path(args.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.use_gpu and args.workers != 1:
        print("--use_gpu: forcing --workers 1 (one CUDA context at a time)")
        args.workers = 1
    sizes = [int(x) for x in floats(args.sizes)] if args.sizes else None
    if args.protocol == "qs":
        jobs = build_qs_jobs(out_dir, sizes=sizes,
                             labels=args.qs_labels.split(",") if args.qs_labels else None,
                             seeds=[int(x) for x in args.qs_seeds.split(",")] if args.qs_seeds else None,
                             gpu=args.use_gpu)
    else:
        jobs = build_jobs(out_dir, sizes=sizes, q0_values=floats(args.q0), r_values=floats(args.rates),
                          gpu=args.use_gpu, seed_check=not args.no_seed_check)
    shard_info = None
    if args.shard:
        k, n = (int(x) for x in args.shard.split("/"))
        weights = floats(args.shard_weights) or [1.0] * n
        if len(weights) != n or not 1 <= k <= n:
            raise SystemExit("--shard k/n needs 1 <= k <= n and n weights")
        n_all = len(jobs)
        jobs = shard_jobs(jobs, k, n, weights)
        shard_info = dict(k=k, n=n, weights=weights, n_jobs_all=n_all)
        print(f"shard {k}/{n} (weights {weights}): {len(jobs)} of {n_all} jobs")
    if args.skip_tags:
        skip = {t.strip() for t in args.skip_tags.split(",") if t.strip()}
        n_before = len(jobs)
        jobs = [j for j in jobs if j["tag"] not in skip]
        print(f"--skip_tags: skipping {n_before - len(jobs)} job(s) done elsewhere")
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
    plan_name = "campaign_plan.json" if not shard_info else f"campaign_plan_shard{shard_info['k']}of{shard_info['n']}.json"
    (out_dir / plan_name).write_text(json.dumps(dict(
        geometries=GEOMETRIES, sizes=SIZES, Q0=Q0_VALUES, r=R_VALUES, r_extra_small=R_EXTRA_SMALL,
        small_N_max=SMALL_N_MAX, seeds=SEEDS, seed_check=SEED_CHECK, B_max_factor=B_MAX_FACTOR,
        needle_thickness=NEEDLE_THICKNESS, n_jobs=len(jobs), est_cpu_h=total / 3600,
        cli=vars(args), shard=shard_info,
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
