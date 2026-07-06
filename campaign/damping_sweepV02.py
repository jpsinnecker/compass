#!/usr/bin/env python3
"""
damping_sweepV02.py
===================

Simulation campaign for the central study:
    "How does the damping coefficient control the shape of the
    hysteresis loop and the statistics of avalanches in classical dipolar
    lattices, and how does this depend on the geometry (square / triangular /
    honeycomb)?"

New in V02 (relative to V01):
    --n_workers N : executes N runs in parallel (independent processes;
                    nearly linear scaling with number of cores).
                    Each worker runs in its own working directory
                    to prevent collision of temporary files
                    (hysteresis_loop.csv) recorded by the simulator in
                    the current directory.
    incremental manifest : the manifest.csv is updated after EACH completed
                    run (not just at the end) — an interruption in the
                    middle of the campaign does not lose progress.
    --resume 1    : resumes an interrupted campaign — runs already
                    registered in manifest.csv are skipped.

This script DOES NOT modify compass.py. It imports the module and calls
relax() programmatically, sweeping:

    damping  x  seed  x  geometry

for the field mode 'hysteresis' (Stage 1) and 'static' B=0 (Stage 2,
free relaxation for domain statistics).

Output
------
For each run, saves a CSV with the complete field_log:
    out_dir/hysteresis/<geom>_damp<idx>_seed<seed>.csv
    out_dir/relax/<geom>_damp<idx>_seed<seed>.csv      (columns: t,B,M,S)

And a metadata JSON per run (derived physical parameters, Q, etc.):
    out_dir/meta/<geom>_damp<idx>_seed<seed>.json

A manifest.csv at the top of out_dir summarizes all runs (one row
per run), updated incrementally.

Usage
-----
    python3 damping_sweepV02.py --out_dir /home/jps/sweep_results \
        --n_seeds 5 --grid_n 30 --n_workers 8

    --quick_test 1   runs a reduced version (fewer dampings/seeds,
                      short t_sim) just to validate that the pipeline
                      works end-to-end before the full campaign.
"""

import os
import sys
import json
import time
import argparse
import tempfile
import itertools
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compass as cs


# ────────────────────────────────────────────────────────────────────────
# Fixed physical parameters of the study (same "tabletop" needle in all
# runs; only geometry/R/damping/seed vary)
# ────────────────────────────────────────────────────────────────────────
R_DEFAULT          = 0.025      # m
NEEDLE_FRAC        = 0.80       # same CLI default
DT_FACTOR          = 0.05       # same CLI default
CUTOFF_DEFAULT     = 3.5        # same CLI default (in units of R)
NOISE_DEFAULT      = 1.5        # rad, initial noise in angles (default make_grid)


def derive_needle_geometry(R, needle_frac=0.80, needle_thickness=0.0004):
    needle_len   = needle_frac * 2.0 * R
    needle_width = needle_len * 0.22
    thickness    = needle_thickness
    return needle_len, needle_width, thickness


def compute_physical_params(R, needle_frac=0.80, needle_thickness=0.0004,
                            steel_density=7850.0, steel_Bsat=None,
                            moment=None, inertia=None):
    """Replicates exactly what main() does to derive inertia/moment
    from geometry, so that the study is physically consistent
    among the 3 geometries (same real needle, same steel density)."""
    needle_len, needle_width, thickness = derive_needle_geometry(R, needle_frac, needle_thickness)
    if inertia is None:
        inertia_val = cs.compute_inertia_from_geometry(needle_len, needle_width, thickness, density=steel_density)
    else:
        inertia_val = inertia

    if moment is None:
        ms_used = steel_Bsat / (4.0 * np.pi * 1e-7) if steel_Bsat is not None else cs.STEEL_MS_SATURATION_DEFAULT
        moment_val = cs.compute_moment_from_geometry(needle_len, needle_width, thickness, Ms=ms_used)
    else:
        moment_val = moment

    return inertia_val, moment_val, needle_len, needle_width, thickness


def compute_Q(inertia, moment, damping, r_nn):
    """Replicates the internal calculation of relax() for Q = omega0*I/b,
    using the dipolar B_ref between nearest neighbors as the scale field."""
    B_ref  = cs.MU0_OVER_4PI * 2.0 * moment / r_nn**3
    omega0 = np.sqrt(moment * B_ref / inertia)
    Q      = omega0 * inertia / damping if damping > 0 else np.inf
    return Q, omega0, B_ref


def make_damping_grid(inertia, moment, r_nn, n_points=8, Q_min=0.05, Q_max=15.0):
    """Constructs the damping grid log-spaced in Q, then converts to
    real `damping` values via the relation Q = omega0*I/b."""
    B_ref  = cs.MU0_OVER_4PI * 2.0 * moment / r_nn**3
    omega0 = np.sqrt(moment * B_ref / inertia)
    Q_vals = np.geomspace(Q_max, Q_min, n_points)   # from most underdamped to most overdamped
    damping_vals = omega0 * inertia / Q_vals
    return damping_vals, Q_vals, omega0


def run_one_hysteresis(geometry, R, damping, seed, t_sim, B_max, n_grid,
                        pbc=False, dt_factor=0.05, use_gpu=False,
                        show_progress=False, damping_noise=0.0,
                        needle_frac=0.80, needle_thickness=0.0004,
                        steel_density=7850.0, steel_Bsat=None,
                        moment=None, inertia=None, noise=1.5):
    np.random.seed(seed)
    inertia_val, moment_val, needle_len, needle_width, thickness = compute_physical_params(
        R, needle_frac=needle_frac, needle_thickness=needle_thickness,
        steel_density=steel_density, steel_Bsat=steel_Bsat,
        moment=moment, inertia=inertia)

    xs, ys, thetas, r_nn, Lx, Ly = cs.make_grid(
        N=n_grid, M=n_grid, geometry=geometry, noise=noise, R=R)

    Q, omega0, B_ref = compute_Q(inertia_val, moment_val, damping, r_nn)

    phi_ext_deg = 0.0  # fixed field direction; symmetry axis chosen per geometry does not matter for M_proj

    t0 = time.perf_counter()
    theta_f, omega_f, hist, n_frames, dt, stop_reason, field_log = cs.relax(
        thetas, xs, ys,
        t_sim=t_sim, dt_factor=dt_factor,
        inertia=inertia_val, damping=damping,
        damping_noise=damping_noise,
        cutoff=CUTOFF_DEFAULT, ext_field=(B_max, 0.0),
        moment=moment_val, field_mode='hysteresis',
        pbc=pbc, Lx=Lx if pbc else None, Ly=Ly if pbc else None,
        use_gpu=use_gpu, show_progress=show_progress,
        make_images=False,
    )
    wall_s = time.perf_counter() - t0

    meta = dict(
        geometry=geometry, R=R, n_grid=n_grid, seed=seed,
        damping=damping, damping_noise=damping_noise, Q=Q, omega0=omega0, B_ref=B_ref, B_max=B_max,
        inertia=inertia_val, moment=moment_val, r_nn=r_nn,
        needle_len=needle_len, needle_width=needle_width, thickness=thickness,
        t_sim=t_sim, dt=dt, dt_factor=dt_factor, n_steps=len(field_log),
        pbc=pbc, stop_reason=stop_reason, wall_seconds=wall_s,
        field_mode='hysteresis',
    )
    return field_log, meta


def run_one_relax(geometry, R, damping, seed, t_sim, n_grid,
                   pbc=False, dt_factor=0.05, use_gpu=False,
                   show_progress=False, damping_noise=0.0,
                   needle_frac=0.80, needle_thickness=0.0004,
                   steel_density=7850.0, steel_Bsat=None,
                   moment=None, inertia=None, noise=1.5,
                   domain_tol=15.0):
    np.random.seed(seed)
    inertia_val, moment_val, needle_len, needle_width, thickness = compute_physical_params(
        R, needle_frac=needle_frac, needle_thickness=needle_thickness,
        steel_density=steel_density, steel_Bsat=steel_Bsat,
        moment=moment, inertia=inertia)

    xs, ys, thetas, r_nn, Lx, Ly = cs.make_grid(
        N=n_grid, M=n_grid, geometry=geometry, noise=noise, R=R)

    Q, omega0, B_ref = compute_Q(inertia_val, moment_val, damping, r_nn)

    t0 = time.perf_counter()
    theta_f, omega_f, hist, n_frames, dt, stop_reason, field_log = cs.relax(
        thetas, xs, ys,
        t_sim=t_sim, dt_factor=dt_factor,
        inertia=inertia_val, damping=damping,
        damping_noise=damping_noise,
        cutoff=CUTOFF_DEFAULT, ext_field=(0.0, 0.0),
        moment=moment_val, field_mode='static',
        pbc=pbc, Lx=Lx if pbc else None, Ly=Ly if pbc else None,
        use_gpu=use_gpu, show_progress=show_progress,
        make_images=False,
    )
    wall_s = time.perf_counter() - t0

    # domain statistics in final state
    domain_labels, n_domains = cs.label_magnetic_domains(theta_f, tol_deg=domain_tol)
    sizes = np.bincount(domain_labels.ravel())
    sizes = sizes[sizes > 0]

    meta = dict(
        geometry=geometry, R=R, n_grid=n_grid, seed=seed,
        damping=damping, damping_noise=damping_noise, Q=Q, omega0=omega0, B_ref=B_ref,
        inertia=inertia_val, moment=moment_val, r_nn=r_nn,
        needle_len=needle_len, needle_width=needle_width, thickness=thickness,
        t_sim=t_sim, dt=dt, dt_factor=dt_factor, n_steps=len(field_log),
        pbc=pbc, stop_reason=stop_reason, wall_seconds=wall_s,
        field_mode='static_relax',
        n_domains=int(n_domains),
        domain_sizes=sizes.tolist(),
    )
    return field_log, meta


def save_field_log_csv(field_log, path):
    arr = np.array(field_log)  # columns: t, B, M, S
    header = "t,B,M,S"
    np.savetxt(path, arr, delimiter=",", header=header, comments="")


# ────────────────────────────────────────────────────────────────────────
# V02: execution by "jobs" — serial or parallel
# ────────────────────────────────────────────────────────────────────────

def _worker_init(scratch_root):
    """Initializer of each worker process: changes to its own working directory
    so that temporary files written by the simulator in the current directory
    (hysteresis_loop.csv) do not collide between concurrent workers."""
    d = tempfile.mkdtemp(prefix="worker_", dir=scratch_root)
    os.chdir(d)


def run_job(job):
    """Executes ONE run (hysteresis or relaxation) and returns the manifest row.
    Module-level function (picklable) for use with ProcessPoolExecutor.
    Writes the CSV and the meta JSON of the run."""
    tag = job["tag"]
    if job["stage"] == "hysteresis":
        field_log, meta = run_one_hysteresis(
            job["geometry"], job["R"], job["damping"], job["seed"],
            job["t_sim"], job["B_max"], job["n_grid"],
            pbc=job["pbc"], use_gpu=job["use_gpu"], show_progress=False,
            damping_noise=job["damping_noise"],
            needle_frac=job["needle_frac"], needle_thickness=job["needle_thickness"],
            steel_density=job["steel_density"], steel_Bsat=job["steel_Bsat"],
            moment=job["moment"], inertia=job["inertia"], dt_factor=job["dt_factor"],
            noise=job["noise"])
        csv_path = os.path.join(job["hyst_dir"], tag + ".csv")
        meta_path = os.path.join(job["meta_dir"], tag + "_hyst.json")
    else:  # relax
        field_log, meta = run_one_relax(
            job["geometry"], job["R"], job["damping"], job["seed"],
            job["t_sim"], job["n_grid"],
            pbc=job["pbc"], use_gpu=job["use_gpu"], show_progress=False,
            damping_noise=job["damping_noise"],
            needle_frac=job["needle_frac"], needle_thickness=job["needle_thickness"],
            steel_density=job["steel_density"], steel_Bsat=job["steel_Bsat"],
            moment=job["moment"], inertia=job["inertia"], dt_factor=job["dt_factor"],
            noise=job["noise"], domain_tol=job["domain_tol"])
        csv_path = os.path.join(job["relax_dir"], tag + ".csv")
        meta_path = os.path.join(job["meta_dir"], tag + "_relax.json")

    save_field_log_csv(field_log, csv_path)
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)

    row = dict(
        tag=tag, stage=job["stage"], geometry=job["geometry"],
        damp_idx=job["damp_idx"], seed=job["seed"], damping=job["damping"],
        Q=meta["Q"], n_steps=meta["n_steps"],
        wall_seconds=meta["wall_seconds"],
        csv_path=csv_path, meta_path=meta_path,
    )
    return row


MANIFEST_FIELDS = ["tag", "stage", "geometry", "damp_idx", "seed",
                   "damping", "Q", "n_steps", "wall_seconds",
                   "csv_path", "meta_path"]


def _append_manifest_row(manifest_path, row):
    """INCREMENTAL write of the manifest: appends the row immediately
    after the run concludes (with flush), so that an interruption in the
    middle of the campaign does not lose progress of already completed runs."""
    import csv as csv_module
    new_file = not os.path.exists(manifest_path)
    with open(manifest_path, "a", newline="") as f:
        writer = csv_module.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        if new_file:
            writer.writeheader()
        writer.writerow(row)
        f.flush()


def _load_done_set(manifest_path):
    """For --resume: reads the existing manifest and returns the set of
    already completed (stage, tag) runs."""
    import csv as csv_module
    done = set()
    if os.path.exists(manifest_path):
        with open(manifest_path, newline="") as f:
            for row in csv_module.DictReader(f):
                done.add((row["stage"], row["tag"]))
    return done


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out_dir', type=str, required=True)
    ap.add_argument('--grid_n', type=int, default=30,
                     help='N=M of the grid (default 30 -> 900 needles)')
    ap.add_argument('--n_seeds', type=int, default=5)
    ap.add_argument('--n_dampings', type=int, default=8)
    ap.add_argument('--Q_min', type=float, default=0.05)
    ap.add_argument('--Q_max', type=float, default=15.0)
    ap.add_argument('--geometries', type=str, default='square,triangular,honeycomb')
    ap.add_argument('--R', type=float, default=R_DEFAULT)
    ap.add_argument('--B_max_factor', type=float, default=8.0,
                     help='B_max = factor * B_ref (dipolar field between neighbors), '
                          'to ensure saturation at the peak of the cycle')
    ap.add_argument('--t_sim_periods', type=float, default=40.0,
                     help='t_sim in units of T0 = 2*pi/omega0 (natural scale). '
                          'Hysteresis cycle has 5 segments; this is the TOTAL time of the cycle.')
    ap.add_argument('--pbc', type=int, default=0)
    ap.add_argument('--use_gpu', type=int, default=0)
    ap.add_argument('--skip_relax_stage', type=int, default=0,
                     help='If 1, only runs Stage 1 (hysteresis), skips Stage 2 (relaxation/domains)')
    ap.add_argument('--damping_noise', type=float, default=0.0,
                     help='Relative noise amplitude for per-needle damping')
    ap.add_argument('--needle_frac', type=float, default=0.80,
                     help='Needle length fraction (0.0 to 0.8)')
    ap.add_argument('--needle_thickness', type=float, default=0.0004,
                     help='Thickness of the steel sheet [m]')
    ap.add_argument('--steel_density', type=float, default=7850.0,
                     help='Density of the steel [kg/m^3]')
    ap.add_argument('--steel_Bsat', type=float, default=None,
                     help='Saturation flux density of the steel [T]')
    ap.add_argument('--moment', type=float, default=None,
                     help='Magnetic moment overriding calculation')
    ap.add_argument('--inertia', type=float, default=None,
                     help='Moment of inertia overriding calculation')
    ap.add_argument('--dt_factor', type=float, default=0.05,
                     help='Fraction of natural period T0 used as time step dt')
    ap.add_argument('--noise', type=float, default=1.5,
                     help='Initial orientation noise amplitude [rad]')
    ap.add_argument('--domain_tol', type=float, default=15.0,
                     help='Angular tolerance [degrees] for domain grouping')
    ap.add_argument('--n_workers', type=int, default=1,
                     help='Number of runs executed in parallel (processes). '
                          '1 = serial (default). Almost linear scaling with the number of cores. '
                          'Avoid combining >1 with --use_gpu 1 (GPU contention).')
    ap.add_argument('--resume', type=int, default=0,
                     help='If 1, resumes interrupted campaign: runs already present '
                          'in manifest.csv of out_dir are skipped.')
    ap.add_argument('--quick_test', type=int, default=0,
                     help='If 1, overrides n_seeds=1, n_dampings=3, t_sim_periods=8, '
                          'grid_n=12 -- just to validate the pipeline end-to-end')
    args = ap.parse_args()

    if args.quick_test:
        args.n_seeds = 1
        args.n_dampings = 3
        args.t_sim_periods = 8.0
        args.grid_n = 12
        print(">>> QUICK_TEST MODE: reduced parameters for pipeline validation <<<")

    geometries = args.geometries.split(',')

    out_dir = args.out_dir
    hyst_dir   = os.path.join(out_dir, 'hysteresis')
    relax_dir  = os.path.join(out_dir, 'relax')
    meta_dir   = os.path.join(out_dir, 'meta')
    for d in (hyst_dir, relax_dir, meta_dir):
        os.makedirs(d, exist_ok=True)

    manifest_path = os.path.join(out_dir, 'manifest.csv')

    # ── damping grid: calculated once with the reference geometry
    #    (square), and reapplied to other geometries. This means that
    #    the SAME set of physical damping values [N.m.s/rad] is used
    #    in all geometries -- what changes between geometries is r_nn (and
    #    therefore effective Q varies slightly between geometries, which will
    #    be recorded in the metadata of each run and used in the analysis, not hidden).
    inertia_ref, moment_ref, _, _, _ = compute_physical_params(
        args.R, needle_frac=args.needle_frac, needle_thickness=args.needle_thickness,
        steel_density=args.steel_density, steel_Bsat=args.steel_Bsat,
        moment=args.moment, inertia=args.inertia)
    xs_ref, ys_ref, _, r_nn_ref, _, _ = cs.make_grid(
        N=args.grid_n, M=args.grid_n, geometry='square', noise=0.0, R=args.R)
    damping_vals, Q_nominal, omega0_ref = make_damping_grid(
        inertia_ref, moment_ref, r_nn_ref,
        n_points=args.n_dampings, Q_min=args.Q_min, Q_max=args.Q_max)

    B_ref_square = cs.MU0_OVER_4PI * 2.0 * moment_ref / r_nn_ref**3
    B_max = args.B_max_factor * B_ref_square
    T0_ref = 2.0 * np.pi / omega0_ref
    t_sim = args.t_sim_periods * T0_ref

    print(f"Damping grid (nominal Q, square geometry as reference):")
    for d, q in zip(damping_vals, Q_nominal):
        print(f"    damping={d:.4e}  N.m.s/rad   ->  Q_nominal={q:.3f}")
    print(f"B_max = {B_max:.5f} T  ({args.B_max_factor}x B_ref)")
    print(f"t_sim = {t_sim:.4f} s  ({args.t_sim_periods} periods T0={T0_ref:.5f}s)")
    print()

    total_runs = len(geometries) * len(damping_vals) * args.n_seeds
    t_campaign_start = time.perf_counter()

    # ── V02: builds the complete list of jobs (Stages 1 and 2) ──────────────
    common = dict(R=args.R, t_sim=t_sim, B_max=B_max, n_grid=args.grid_n,
                  pbc=bool(args.pbc), use_gpu=bool(args.use_gpu),
                  damping_noise=args.damping_noise,
                  needle_frac=args.needle_frac, needle_thickness=args.needle_thickness,
                  steel_density=args.steel_density, steel_Bsat=args.steel_Bsat,
                  moment=args.moment, inertia=args.inertia, dt_factor=args.dt_factor,
                  noise=args.noise, domain_tol=args.domain_tol,
                  hyst_dir=hyst_dir, relax_dir=relax_dir, meta_dir=meta_dir)
    jobs = []
    stages = ['hysteresis'] + ([] if args.skip_relax_stage else ['relax'])
    for stage in stages:
        for geometry in geometries:
            for damp_idx, damping in enumerate(damping_vals):
                for seed in range(args.n_seeds):
                    tag = f"{geometry}_damp{damp_idx:02d}_seed{seed:02d}"
                    jobs.append(dict(common, stage=stage, tag=tag,
                                     geometry=geometry, damp_idx=damp_idx,
                                     seed=seed, damping=damping))

    # ── V02: resume — skips jobs already registered in manifest ──────────────
    if args.resume:
        done = _load_done_set(manifest_path)
        n_before = len(jobs)
        jobs = [j for j in jobs if (j['stage'], j['tag']) not in done]
        print(f"Resume: {n_before - len(jobs)} runs already completed "
              f"(skipped); {len(jobs)} remaining.")
    elif os.path.exists(manifest_path):
        # new campaign over old out_dir: starts the manifest from scratch
        os.remove(manifest_path)

    n_total = len(jobs)
    if n_total == 0:
        print("Nothing to do (all runs already completed).")
        return

    # ── V02: execution — serial or parallel — with incremental manifest ───
    scratch_root = os.path.join(out_dir, '_scratch')
    os.makedirs(scratch_root, exist_ok=True)

    if args.n_workers <= 1:
        # serial: runs in current working directory (V01 behavior)
        for i, job in enumerate(jobs, 1):
            print(f"[{i}/{n_total}] {job['stage'].upper():10s} {job['tag']}  "
                  f"(damping={job['damping']:.3e})  ...", end=' ', flush=True)
            row = run_job(job)
            _append_manifest_row(manifest_path, row)
            print(f"ok  ({row['wall_seconds']:.1f}s, Q={row['Q']:.3f})")
    else:
        print(f"Running with {args.n_workers} workers in parallel...")
        n_done = 0
        with ProcessPoolExecutor(max_workers=args.n_workers,
                                 initializer=_worker_init,
                                 initargs=(scratch_root,)) as ex:
            futures = {ex.submit(run_job, job): job for job in jobs}
            for fut in as_completed(futures):
                job = futures[fut]
                row = fut.result()   # propagates exceptions from worker
                _append_manifest_row(manifest_path, row)
                n_done += 1
                print(f"[{n_done}/{n_total}] {job['stage'].upper():10s} "
                      f"{job['tag']}  ok  ({row['wall_seconds']:.1f}s, "
                      f"Q={row['Q']:.3f})", flush=True)

    total_wall = time.perf_counter() - t_campaign_start
    print()
    print(f"Campaign completed in {total_wall/60.0:.1f} min.")
    print(f"Manifest: {manifest_path}")


if __name__ == '__main__':
    main()
