#!/usr/bin/env python3
"""
sweep_damping_hysteresis.py
============================

NOTE (added after the fact): this script is non-functional outside
--dry_run -- run_single() references undefined module globals (NameError),
miscalls compute_t_sim() (a positional arg silently binds to the wrong
parameter), and passes a damp_idx kwarg run_single() doesn't accept
(TypeError). See docs/AUDIT.md bug B3. Use damping_sweepV03.py instead
(see USER_GUIDE.md Sec 9). Left as-is, not patched, per docs/AUDIT.md's
own P0 recommendation.

Sweep campaign: damping (Q) x geometry x seed, in --field_mode hysteresis.

Implemented script:
  - 8 log-spaced damping values covering Q ~ 15 (underdamped) to Q ~ 0.05
    (overdamped), calculated from the real physical system (omega0, I) —
    not arbitrary values.
  - 3 geometries: square, triangular, honeycomb
  - 5 initial noise seeds per (damping, geometry) combination
  - Each run executes a full hysteresis cycle (5 segments:
    0->+Bmax->0->-Bmax->0->+Bmax) and saves the field_log (t, B, M, S)
    to CSV, plus a metadata JSON for the run.

Total: 8 x 3 x 5 = 120 runs.

Usage:
    python3 sweep_damping_hysteresis.py --out_dir /home/jps/sweep_results
    python3 sweep_damping_hysteresis.py --out_dir ... --dry_run   # only list the plan
    python3 sweep_damping_hysteresis.py --out_dir ... --resume    # skip already existing CSVs
    python3 sweep_damping_hysteresis.py --out_dir ... --gpu 1     # use GPU if available

Each run is independent -> can be interrupted (Ctrl+C) and resumed
later with --resume without losing already completed work.
"""

import argparse
import contextlib
import csv
import io
import json
import os
import shutil
import sys
import time
import traceback

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compass as cs

from sim_config import load_config

_CFG = load_config()
_PHYS = _CFG.physics.sweep_damping_hysteresis
_NUM = _CFG.numerics.sweep_damping_hysteresis
_RUN = _CFG.run.sweep_damping_hysteresis


# ─────────────────────────────────────────────────────────────────────────
# Reference physical parameters — derived from needle geometry using the
# same functions as compass.py main() and damping_sweep.py V01/V02.
# Using cs.MOMENT_DEFAULT/INERTIA_DEFAULT is WRONG here because those are
# fallback values, not the geometry-derived ones for R=0.025 m needle.
# ─────────────────────────────────────────────────────────────────────────
R_DEFAULT       = _PHYS.R_default      # m
NEEDLE_FRAC     = _PHYS.needle_frac    # same as compass.py default
NEEDLE_THICKNESS = _PHYS.needle_thickness  # m
STEEL_DENSITY   = _PHYS.steel_density  # kg/m³
MU0_OVER_4PI    = cs.MU0_OVER_4PI

N_GRID          = _NUM.N_grid          # N x M -> 900 needles
M_GRID          = _NUM.M_grid
NOISE_INIT      = _PHYS.noise_init     # rad, initial noise in angles
SEEDS           = _RUN.seeds
GEOMETRIES      = _RUN.geometries
N_DAMPING       = _NUM.n_damping
Q_MIN, Q_MAX    = _NUM.Q_min, _NUM.Q_max  # overdamped -> underdamped

# Interaction cutoff in units of r_nn (nearest-neighbour distance = 2R).
# 3.5 matches compass.py CLI default and damping_sweep.py V01/V02.
# (The old value of 8.0 was metres, not r_nn units, and is incorrect.)
CUTOFF          = _NUM.cutoff

# B_max scale factor: B_max = B_MAX_FACTOR * B_ref (dipolar field between
# nearest neighbours). 8.0 matches damping_sweep.py V01/V02 default.
B_MAX_FACTOR    = _PHYS.B_max_factor

DT_FACTOR       = _NUM.dt_factor       # fraction of T0 per step
T_SIM_IN_T0     = _NUM.t_sim_in_T0     # total cycle duration in units of T0


def _compute_derived_params(R=R_DEFAULT):
    """Derive needle geometry and moment/inertia from physical parameters,
    consistent with compass.py main() and damping_sweep.py V01/V02."""
    needle_len   = NEEDLE_FRAC * 2.0 * R
    needle_width = needle_len * _CFG.physics.needle_geometry.default_width_to_length_ratio
    thickness    = NEEDLE_THICKNESS
    inertia = cs.compute_inertia_from_geometry(
        needle_len, needle_width, thickness, density=STEEL_DENSITY,
        pivot_radius=_PHYS.pivot_radius, pivot_thickness=_PHYS.pivot_thickness)
    moment  = cs.compute_moment_from_geometry(
        needle_len, needle_width, thickness, Ms=cs.STEEL_MS_SATURATION_DEFAULT,
        pivot_radius=_PHYS.pivot_radius)
    return moment, inertia, needle_len, needle_width, thickness


def compute_omega0(R=R_DEFAULT):
    """Replicates the omega0 calculation done inside relax() (B_ref from the
    dipolar field between nearest neighbours), to define the damping range
    and t_sim consistently with the real physical system simulated."""
    moment, inertia, _, _, _ = _compute_derived_params(R)
    r_nn = 2.0 * R  # all geometries have r_nn = 2R (nearest-neighbour distance)
    B_ref = MU0_OVER_4PI * 2.0 * moment / r_nn**3
    omega0 = np.sqrt(moment * B_ref / inertia)
    return omega0, moment, inertia


def compute_damping_values(R=R_DEFAULT, n=N_DAMPING, q_min=Q_MIN, q_max=Q_MAX):
    """8 log-spaced damping values covering Q from q_max (underdamped)
    to q_min (overdamped). Returns a list of dicts {Q, damping}."""
    omega0, moment, inertia = compute_omega0(R)
    Q_targets = np.logspace(np.log10(q_max), np.log10(q_min), n)
    b_values = omega0 * inertia / Q_targets
    return [{"Q": float(Q), "damping": float(b)} for Q, b in zip(Q_targets, b_values)]


def compute_t_sim(R=R_DEFAULT):
    """t_sim in seconds: fix cycle duration in units of the natural period T0."""
    omega0, _, _ = compute_omega0(R)
    T0 = 2.0 * np.pi / omega0
    return T_SIM_IN_T0 * T0


def run_single(geometry, damping, Q, seed, out_dir, R=R_DEFAULT,
                use_gpu=False, show_progress=False):
    """Runs a hysteresis simulation and saves the field_log + metadata.
    Returns the saved CSV path and wall time."""

    tag = f"{geometry}_Q{Q:08.4f}_seed{seed:02d}"
    csv_path  = os.path.join(out_dir, f"{tag}.csv")
    meta_path = os.path.join(out_dir, f"{tag}.json")

    np.random.seed(seed)
    xs, ys, thetas, nn_dist, Lx, Ly = cs.make_grid(
        N=N_GRID, M=M_GRID, geometry=geometry, noise=NOISE_INIT, R=R
    )

    t_sim = compute_t_sim(damping, R=R)

    # relax() writes "hysteresis_loop.csv"/".png" in the cwd as a hardcoded side
    # effect (in addition to returning field_log, which is what we actually use).
    # We isolate each run in its own scratch dir to avoid overwriting/accumulating
    # garbage, and discard this scratch at the end -- the "official" CSV of this
    # run is what we write below, in out_dir.
    scratch_dir = os.path.join(out_dir, "_scratch", tag)
    os.makedirs(scratch_dir, exist_ok=True)
    cwd_orig = os.getcwd()

    t_wall_start = time.perf_counter()
    try:
        os.chdir(scratch_dir)
        # show_progress=False still leaves informative _print() calls from
        # relax() active (there is no verbosity flag in compass.py) --
        # we suppress stdout during the call to avoid polluting the log of
        # the 120-run campaign.
        _buf = io.StringIO()
        with contextlib.redirect_stdout(_buf):
            theta_cur, omega_cur, hist, n_frames, dt, stop_reason, field_log = cs.relax(
                thetas, xs, ys,
                t_sim=t_sim,
                dt_factor=DT_FACTOR,
                inertia=INERTIA,
                damping=damping,
                cutoff=CUTOFF,
                ext_field=(B_MAX, 0.0),
                moment=MOMENT,
                field_mode="hysteresis",
                frame_dir=None,
                pbc=False,
                use_gpu=use_gpu,
                show_progress=show_progress,
            )
    finally:
        os.chdir(cwd_orig)
        shutil.rmtree(scratch_dir, ignore_errors=True)
    t_wall = time.perf_counter() - t_wall_start

    # ── saves the field_log (t, B_proj, M_proj, S) in CSV ────────────────────
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["t", "B_proj", "M_proj", "S"])
        writer.writerows(field_log)

    # ── run metadata (essential to reproduce/audit later) ────
    meta = {
        "geometry": geometry,
        "Q_target": Q,
        "damping": damping,
        "seed": seed,
        "R": R,
        "N_grid": N_GRID,
        "M_grid": M_GRID,
        "noise_init": NOISE_INIT,
        "cutoff": CUTOFF,
        "B_max": B_MAX,
        "dt_factor": DT_FACTOR,
        "dt_actual": dt,
        "t_sim": t_sim,
        "n_steps_logged": len(field_log),
        "stop_reason": stop_reason,
        "wall_time_seconds": t_wall,
        "moment": MOMENT,
        "inertia": INERTIA,
        "used_gpu": bool(use_gpu),
    }
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)

    return csv_path, t_wall


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out_dir", required=True, help="Output directory for CSVs/JSONs")
    ap.add_argument("--dry_run", action="store_true", help="Only lists the run plan, does not execute")
    ap.add_argument("--resume", action="store_true", help="Skips runs whose CSV already exists")
    ap.add_argument("--gpu", type=int, default=_RUN.gpu, help="1 to use GPU (CuPy) if available")
    ap.add_argument("--geometries", nargs="+", default=GEOMETRIES, choices=GEOMETRIES)
    ap.add_argument("--seeds", nargs="+", type=int, default=SEEDS)
    ap.add_argument("--n_damping", type=int, default=N_DAMPING)
    ap.add_argument("--only_index", type=int, default=None,
                     help="Runs only the run of index N of the plan (debug/manual parallelization)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    damping_values = compute_damping_values(n=args.n_damping)

    plan = []
    for geometry in args.geometries:
        for damp_idx, dv in enumerate(damping_values):
            for seed in args.seeds:
                plan.append({"geometry": geometry, "Q": dv["Q"],
                             "damp_idx": damp_idx,
                             "damping": dv["damping"], "seed": seed})

    print(f"Plan: {len(plan)} runs "
          f"({len(args.geometries)} geometries x {len(damping_values)} dampings x {len(args.seeds)} seeds)")
    print(f"Q range: {damping_values[-1]['Q']:.3f} (overdamped) "
          f"to {damping_values[0]['Q']:.3f} (underdamped)")
    print(f"t_sim varies by damping (proportional to {T_SIM_IN_T0} x T0)")
    print()

    if args.dry_run:
        for i, p in enumerate(plan):
            print(f"  [{i:3d}] geometry={p['geometry']:<10s} Q={p['Q']:8.4f} "
                  f"damping={p['damping']:.3e} seed={p['seed']}")
        return

    if args.only_index is not None:
        plan = [plan[args.only_index]]

    n_done, n_skipped, n_failed = 0, 0, 0
    t_campaign_start = time.perf_counter()

    for i, p in enumerate(plan):
        tag = f"{p['geometry']}_Q{p['Q']:08.4f}_seed{p['seed']:02d}"
        csv_path = os.path.join(args.out_dir, f"{tag}.csv")

        if args.resume and os.path.exists(csv_path):
            n_skipped += 1
            print(f"[{i+1}/{len(plan)}] {tag}  -- already exists, skipping (--resume)")
            continue

        print(f"[{i+1}/{len(plan)}] {tag}  -- running...", flush=True)
        try:
            _, t_wall = run_single(
                geometry=p["geometry"], damping=p["damping"], Q=p["Q"],
                seed=p["seed"], out_dir=args.out_dir, damp_idx=p["damp_idx"],
                use_gpu=bool(args.gpu), show_progress=False,
            )
            n_done += 1
            print(f"           ok  ({t_wall:.1f} s)")
        except KeyboardInterrupt:
            print("\nInterrupted by user (Ctrl+C). "
                  "Use --resume to continue from where you stopped.")
            sys.exit(130)
        except Exception as e:
            n_failed += 1
            print(f"           FAILED: {e}")
            traceback.print_exc()
            # logs the failure to not silence problems
            fail_path = os.path.join(args.out_dir, f"{tag}.FAILED.txt")
            with open(fail_path, "w") as f:
                f.write(f"{type(e).__name__}: {e}\n\n")
                traceback.print_exc(file=f)

    t_campaign = time.perf_counter() - t_campaign_start
    print()
    print(f"Campaign completed in {t_campaign/60.0:.1f} min: "
          f"{n_done} ok, {n_skipped} skipped, {n_failed} failed.")


if __name__ == "__main__":
    main()
