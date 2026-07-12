#!/usr/bin/env python3
"""
damping_sweep.py
=================

Simulation campaign for the central study:
    "How does the damping coefficient control the shape of the
    hysteresis loop and the statistics of avalanches in classical dipolar
    lattices, and how does this depend on the geometry (square / triangular /
    honeycomb)?"

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
per run), to facilitate loading in pandas during analysis.

Usage
-----
    python3 damping_sweep.py --out_dir /home/jps/sweep_results \
        --n_seeds 5 --grid_n 30 --quick_test 1

    --quick_test 1   runs a reduced version (fewer dampings/seeds,
                      short t_sim) just to validate that the pipeline
                      works end-to-end before the full campaign.
    --quick_test 0   runs the full campaign according to the script.
"""

import os
import sys
import json
import time
import argparse
import itertools
import contextlib
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compass as cs

from sim_config import load_config

_CFG = load_config()
_PHYS = _CFG.physics.damping_sweep
_NUM = _CFG.numerics.damping_sweep
_RUN = _CFG.run.damping_sweep


# ────────────────────────────────────────────────────────────────────────
# Fixed physical parameters of the study (same "tabletop" needle in all
# runs; only geometry/R/damping/seed vary)
# ────────────────────────────────────────────────────────────────────────
# Real compass-needle geometry (small tabletop apparatus):
#   R = 7.5 mm, needle_frac = 2/3, thickness = 0.26 mm,
#   pivot pin: r=1 mm, h=2 mm, brass (8500 kg/m³)
# For the LAW3M 2026 campaign (large demonstration needle):
#   override with --R 0.025 --needle_frac 0.80 --needle_thickness 0.4e-3
#                  --pivot_radius 0 --pivot_thickness 0
R_DEFAULT          = _PHYS.R_default      # m  — real small needle
NEEDLE_FRAC        = _PHYS.needle_frac    # same CLI default
DT_FACTOR          = _NUM.dt_factor       # same CLI default
CUTOFF_DEFAULT     = _PHYS.cutoff_default  # in units of r_nn = 2R
NOISE_DEFAULT      = _PHYS.noise          # rad, initial orientation noise


def _bar(done, total, label="", extra="", width=36):
    """Print a one-line overwriting progress bar.
    Call with done=0 before the loop to show an empty bar, and with
    done==total after the last job to finalise (adds newline).

    Example output (on a single overwritten line):
      HYSTERESIS  [████████████░░░░░░░░░░░░░░░░░░░░░░]  12/80  15.0%  Q=3.46
    """
    frac   = done / total if total > 0 else 0.0
    filled = int(width * frac)
    bar    = "█" * filled + "░" * (width - filled)
    pct    = frac * 100.0
    suffix = f"  {extra}" if extra else ""
    line   = f"  {label:10s} [{bar}] {done:>4d}/{total}  {pct:5.1f}%{suffix}"
    end    = "\n" if done == total else "\r"
    print(line, end=end, flush=True)


def derive_needle_geometry(R, needle_frac=2.0/3.0, needle_thickness=0.26e-3):
    needle_len   = needle_frac * 2.0 * R
    needle_width = needle_len * _CFG.physics.needle_geometry.default_width_to_length_ratio
    thickness    = needle_thickness
    return needle_len, needle_width, thickness


def compute_physical_params(R, needle_frac=2.0/3.0, needle_thickness=0.26e-3,
                            steel_density=7850.0, steel_Bsat=None,
                            moment=None, inertia=None,
                            pivot_radius=1.0e-3, pivot_thickness=2.0e-3,
                            pivot_density=8500.0, pivot_mass=None):
    """Replicates exactly what main() does to derive inertia/moment
    from geometry, so that the study is physically consistent
    among the 3 geometries (same real needle, same steel density).

    Pivot parameters (pivot_radius, pivot_thickness, pivot_density, pivot_mass)
    are forwarded to compute_inertia_from_geometry() to match compass.py V78's
    full inertia calculation, which includes the physical pivot pin.
    All default to zero (no pivot), consistent with compass.py defaults.
    """
    needle_len, needle_width, thickness = derive_needle_geometry(R, needle_frac, needle_thickness)
    if inertia is None:
        inertia_val = cs.compute_inertia_from_geometry(
            needle_len, needle_width, thickness, density=steel_density,
            pivot_radius=pivot_radius, pivot_thickness=pivot_thickness,
            pivot_density=pivot_density, pivot_mass=pivot_mass)
    else:
        inertia_val = inertia

    if moment is None:
        ms_used = steel_Bsat / (4.0 * np.pi * 1e-7) if steel_Bsat is not None else cs.STEEL_MS_SATURATION_DEFAULT
        moment_val = cs.compute_moment_from_geometry(needle_len, needle_width, thickness, Ms=ms_used, pivot_radius=pivot_radius)
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
                        moment=None, inertia=None, noise=1.5,
                        pivot_radius=0.0, pivot_thickness=0.0,
                        pivot_density=8500.0, pivot_mass=None):
    np.random.seed(seed)
    inertia_val, moment_val, needle_len, needle_width, thickness = compute_physical_params(
        R, needle_frac=needle_frac, needle_thickness=needle_thickness,
        steel_density=steel_density, steel_Bsat=steel_Bsat,
        moment=moment, inertia=inertia,
        pivot_radius=pivot_radius, pivot_thickness=pivot_thickness,
        pivot_density=pivot_density, pivot_mass=pivot_mass)

    xs, ys, thetas, r_nn, Lx, Ly = cs.make_grid(
        N=n_grid, M=n_grid, geometry=geometry, noise=noise, R=R)

    Q, omega0, B_ref = compute_Q(inertia_val, moment_val, damping, r_nn)

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
        pivot_radius=pivot_radius, pivot_thickness=pivot_thickness,
        pivot_density=pivot_density, pivot_mass=pivot_mass,
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
                   domain_tol=15.0,
                   pivot_radius=0.0, pivot_thickness=0.0,
                   pivot_density=8500.0, pivot_mass=None):
    np.random.seed(seed)
    inertia_val, moment_val, needle_len, needle_width, thickness = compute_physical_params(
        R, needle_frac=needle_frac, needle_thickness=needle_thickness,
        steel_density=steel_density, steel_Bsat=steel_Bsat,
        moment=moment, inertia=inertia,
        pivot_radius=pivot_radius, pivot_thickness=pivot_thickness,
        pivot_density=pivot_density, pivot_mass=pivot_mass)

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
        pivot_radius=pivot_radius, pivot_thickness=pivot_thickness,
        pivot_density=pivot_density, pivot_mass=pivot_mass,
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


def save_forc_csv(field_log, path):
    """Saves the FORC field log (7 columns) to CSV.
    Columns: t, B, M, S, B_r, sweep_dir, curve_idx
    Only 'up'-sweep rows are written (ascending branches only)."""
    import csv as csv_mod
    with open(path, 'w', newline='') as f:
        w = csv_mod.writer(f)
        w.writerow(['t', 'B', 'M', 'S', 'B_r', 'sweep_dir', 'curve_idx'])
        for row in field_log:
            if len(row) >= 7 and row[5] == 'up':
                w.writerow(row)


def run_one_forc(geometry, R, damping, seed, B_max, n_grid,
                  pbc=False, dt_factor=0.05, use_gpu=False,
                  show_progress=False, damping_noise=0.0,
                  needle_frac=2.0/3.0, needle_thickness=0.26e-3,
                  steel_density=7850.0, steel_Bsat=None,
                  moment=None, inertia=None, noise=1.5,
                  pivot_radius=1.0e-3, pivot_thickness=2.0e-3,
                  pivot_density=8500.0, pivot_mass=None,
                  forc_n_curves=20, forc_Br_min=None,
                  forc_rate=None, forc_t_sat=None,
                  forc_t_ramp_down=None, forc_t_ramp_up=None):
    """Runs one complete FORC measurement (all reversal curves) for a single
    (geometry, damping, seed) point, using compass.py field_mode='forc'.
    Returns (field_log, meta) where field_log rows are
    [t, B, M, S, B_r, sweep_dir, curve_idx]."""
    np.random.seed(seed)
    inertia_val, moment_val, needle_len, needle_width, thickness = compute_physical_params(
        R, needle_frac=needle_frac, needle_thickness=needle_thickness,
        steel_density=steel_density, steel_Bsat=steel_Bsat,
        moment=moment, inertia=inertia,
        pivot_radius=pivot_radius, pivot_thickness=pivot_thickness,
        pivot_density=pivot_density, pivot_mass=pivot_mass)

    xs, ys, thetas, r_nn, Lx, Ly = cs.make_grid(
        N=n_grid, M=n_grid, geometry=geometry, noise=noise, R=R)

    Q, omega0, B_ref = compute_Q(inertia_val, moment_val, damping, r_nn)
    Br_min = forc_Br_min if forc_Br_min is not None else -B_max

    t0 = time.perf_counter()
    theta_f, omega_f, hist, n_frames, dt, stop_reason, field_log = cs.relax(
        thetas, xs, ys,
        t_sim=0.0,              # t_sim is auto-computed by relax() in forc mode
        dt_factor=dt_factor,
        inertia=inertia_val, damping=damping,
        damping_noise=damping_noise,
        cutoff=CUTOFF_DEFAULT, ext_field=(B_max, 0.0),
        moment=moment_val, field_mode='forc',
        pbc=pbc, Lx=Lx if pbc else None, Ly=Ly if pbc else None,
        use_gpu=use_gpu, show_progress=show_progress,
        make_images=False,
        forc_n_curves=forc_n_curves,
        forc_Br_min=Br_min,
        forc_rate=forc_rate,
        forc_t_sat=forc_t_sat,
        forc_t_ramp_down=forc_t_ramp_down,
        forc_t_ramp_up=forc_t_ramp_up,
        forc_sweep='up',        # record only the ascending branch (standard FORC)
    )
    wall_s = time.perf_counter() - t0

    meta = dict(
        geometry=geometry, R=R, n_grid=n_grid, seed=seed,
        damping=damping, damping_noise=damping_noise, Q=Q, omega0=omega0, B_ref=B_ref, B_max=B_max,
        inertia=inertia_val, moment=moment_val, r_nn=r_nn,
        needle_len=needle_len, needle_width=needle_width, thickness=thickness,
        pivot_radius=pivot_radius, pivot_thickness=pivot_thickness,
        pivot_density=pivot_density, pivot_mass=pivot_mass,
        dt=dt, dt_factor=dt_factor, n_steps=len(field_log),
        pbc=pbc, stop_reason=stop_reason, wall_seconds=wall_s,
        field_mode='forc',
        forc_n_curves=forc_n_curves,
        forc_Br_min=Br_min,
    )
    return field_log, meta


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out_dir', type=str, required=True)
    ap.add_argument('--grid_n', type=int, default=_NUM.grid_n,
                     help='N=M of the grid (default 30 -> 900 needles)')
    ap.add_argument('--n_seeds', type=int, default=_NUM.n_seeds)
    ap.add_argument('--n_dampings', type=int, default=_NUM.n_dampings)
    ap.add_argument('--Q_min', type=float, default=_NUM.Q_min)
    ap.add_argument('--Q_max', type=float, default=_NUM.Q_max)
    ap.add_argument('--geometries', type=str, default=_RUN.geometries)
    ap.add_argument('--R', type=float, default=R_DEFAULT)
    ap.add_argument('--B_max_factor', type=float, default=_PHYS.B_max_factor,
                     help='B_max = factor * B_ref (dipolar field between neighbors), '
                          'to ensure saturation at the peak of the cycle')
    ap.add_argument('--t_sim_periods', type=float, default=_NUM.t_sim_periods,
                     help='t_sim in units of T0 = 2*pi/omega0 (natural scale). '
                          'Hysteresis cycle has 5 segments; this is the TOTAL time of the cycle.')
    ap.add_argument('--pbc', type=int, default=_RUN.pbc)
    ap.add_argument('--use_gpu', type=int, default=_RUN.use_gpu)
    ap.add_argument('--damping_noise', type=float, default=_RUN.damping_noise,
                     help='Relative noise amplitude for per-needle damping')
    ap.add_argument('--skip_relax_stage', type=int, default=_RUN.skip_relax_stage,
                     help='If 1, only runs Stage 1 (hysteresis), skips Stage 2 (relaxation/domains)')
    ap.add_argument('--skip_forc_stage', type=int, default=_RUN.skip_forc_stage,
                     help='If 1, skips Stage 3 (FORC measurement)')
    ap.add_argument('--forc_n_curves', type=int, default=_NUM.forc_n_curves,
                     help='Number of FORC reversal curves per run (default 20; quick_test uses 3)')
    ap.add_argument('--forc_Br_min', type=float, default=None,
                     help='Minimum reversal field [T] (default: -B_max, i.e. full range)')
    ap.add_argument('--forc_rate', type=float, default=None,
                     help='Constant field sweep rate for FORC ramps [T/s]. '
                          'If None, fixed forc_t_ramp_down/up times are used.')
    ap.add_argument('--forc_t_sat', type=float, default=None,
                     help='Saturation hold time per FORC cycle [s] (default: compass.py default 0.05 s)')
    ap.add_argument('--forc_t_ramp_down', type=float, default=None,
                     help='Ramp-down time per FORC cycle [s] (default: compass.py default 0.10 s)')
    ap.add_argument('--forc_t_ramp_up', type=float, default=None,
                     help='Ramp-up time per FORC cycle [s] (default: compass.py default 0.20 s)')
    ap.add_argument('--needle_frac', type=float, default=_PHYS.needle_frac,
                     help='Needle length as fraction of 2R. '
                          'Real small needle: 2/3. '
                          'LAW3M large needle: 0.80')
    ap.add_argument('--needle_thickness', type=float, default=_PHYS.needle_thickness,
                     help='Steel blade thickness [m]. '
                          'Real small needle: 0.26e-3. '
                          'LAW3M large needle: 0.4e-3')
    ap.add_argument('--steel_density', type=float, default=_PHYS.steel_density,
                     help='Density of the steel [kg/m^3]')
    ap.add_argument('--steel_Bsat', type=float, default=None,
                     help='Saturation flux density of the steel [T]')
    ap.add_argument('--moment', type=float, default=None,
                     help='Magnetic moment overriding calculation')
    ap.add_argument('--inertia', type=float, default=None,
                     help='Moment of inertia overriding calculation')
    ap.add_argument('--pivot_radius', type=float, default=_PHYS.pivot_radius,
                     help='Pivot pin radius [m] (real needle: 1 mm; '
                          'LAW3M large needle: 0 to ignore pivot)')
    ap.add_argument('--pivot_thickness', type=float, default=_PHYS.pivot_thickness,
                     help='Pivot pin thickness [m] (real needle: 2 mm; '
                          'LAW3M large needle: 0 to ignore pivot)')
    ap.add_argument('--pivot_density', type=float, default=_PHYS.pivot_density,
                     help='Density of the pivot pin material [kg/m^3]')
    ap.add_argument('--pivot_mass', type=float, default=None,
                     help='Override pivot mass directly [kg] (overrides pivot_radius/thickness/density)')
    ap.add_argument('--dt_factor', type=float, default=_NUM.dt_factor,
                     help='Fraction of natural period T0 used as time step dt')
    ap.add_argument('--noise', type=float, default=_PHYS.noise,
                     help='Initial orientation noise amplitude [rad]')
    ap.add_argument('--domain_tol', type=float, default=_NUM.domain_tol,
                     help='Angular tolerance [degrees] for domain grouping')
    ap.add_argument('--quick_test', type=int, default=_RUN.quick_test,
                     help='If 1, overrides n_seeds=1, n_dampings=3, t_sim_periods=8, '
                          'grid_n=12 -- just to validate the pipeline end-to-end')
    args = ap.parse_args()

    if args.quick_test:
        qt = _NUM.quick_test
        args.n_seeds = qt.n_seeds
        args.n_dampings = qt.n_dampings
        args.t_sim_periods = qt.t_sim_periods
        args.grid_n = qt.grid_n
        args.forc_n_curves = qt.forc_n_curves
        print(">>> QUICK_TEST MODE: reduced parameters for pipeline validation <<<")

    geometries = args.geometries.split(',')

    out_dir = args.out_dir
    hyst_dir   = os.path.join(out_dir, 'hysteresis')
    relax_dir  = os.path.join(out_dir, 'relax')
    forc_dir   = os.path.join(out_dir, 'forc')
    meta_dir   = os.path.join(out_dir, 'meta')
    for d in (hyst_dir, relax_dir, forc_dir, meta_dir):
        os.makedirs(d, exist_ok=True)

    manifest_path = os.path.join(out_dir, 'manifest.csv')
    manifest_rows = []

    # ── damping grid: calculated once with the reference geometry
    #    (square), and reapplied to other geometries. This means that
    #    the SAME set of physical damping values [N.m.s/rad] is used
    #    in all geometries -- what changes between geometries is r_nn (and
    #    therefore effective Q varies slightly between geometries, which will
    #    be recorded in the metadata of each run and used in the analysis, not hidden).
    inertia_ref, moment_ref, _, _, _ = compute_physical_params(
        args.R, needle_frac=args.needle_frac, needle_thickness=args.needle_thickness,
        steel_density=args.steel_density, steel_Bsat=args.steel_Bsat,
        moment=args.moment, inertia=args.inertia,
        pivot_radius=args.pivot_radius, pivot_thickness=args.pivot_thickness,
        pivot_density=args.pivot_density, pivot_mass=args.pivot_mass)
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
    run_idx = 0
    t_campaign_start = time.perf_counter()

    # ── Stage 1: hysteresis ──────────────────────────────────────────────
    _bar(0, total_runs, label="HYSTERESIS")
    for geometry in geometries:
        for damp_idx, damping in enumerate(damping_vals):
            for seed in range(args.n_seeds):
                run_idx += 1
                tag = f"{geometry}_damp{damp_idx:02d}_seed{seed:02d}"

                with open(os.devnull, 'w') as _dn, \
                        contextlib.redirect_stdout(_dn):
                    field_log, meta = run_one_hysteresis(
                        geometry, args.R, damping, seed, t_sim, B_max, args.grid_n,
                        pbc=bool(args.pbc), use_gpu=bool(args.use_gpu),
                        damping_noise=args.damping_noise,
                        needle_frac=args.needle_frac, needle_thickness=args.needle_thickness,
                        steel_density=args.steel_density, steel_Bsat=args.steel_Bsat,
                        moment=args.moment, inertia=args.inertia,
                        dt_factor=args.dt_factor, noise=args.noise,
                        pivot_radius=args.pivot_radius, pivot_thickness=args.pivot_thickness,
                        pivot_density=args.pivot_density, pivot_mass=args.pivot_mass)

                csv_path = os.path.join(hyst_dir, tag + '.csv')
                meta_path = os.path.join(meta_dir, tag + '_hyst.json')
                save_field_log_csv(field_log, csv_path)
                with open(meta_path, 'w') as f:
                    json.dump(meta, f, indent=2)

                manifest_rows.append(dict(
                    tag=tag, stage='hysteresis', geometry=geometry,
                    damp_idx=damp_idx, seed=seed, damping=damping,
                    Q=meta['Q'], n_steps=meta['n_steps'],
                    wall_seconds=meta['wall_seconds'],
                    csv_path=csv_path, meta_path=meta_path,
                ))
                _bar(run_idx, total_runs, label="HYSTERESIS",
                     extra=f"{tag}  Q={meta['Q']:.3f}  {meta['wall_seconds']:.1f}s")

    # ── Stage 2: free relaxation / domains ──────────────────────────────
    if not args.skip_relax_stage:
        t_sim_relax = t_sim  # same total time scale, no field
        total_runs_relax = len(geometries) * len(damping_vals) * args.n_seeds
        run_idx2 = 0
        _bar(0, total_runs_relax, label="RELAX")
        for geometry in geometries:
            for damp_idx, damping in enumerate(damping_vals):
                for seed in range(args.n_seeds):
                    run_idx2 += 1
                    tag = f"{geometry}_damp{damp_idx:02d}_seed{seed:02d}"

                    with open(os.devnull, 'w') as _dn, \
                            contextlib.redirect_stdout(_dn):
                        field_log, meta = run_one_relax(
                            geometry, args.R, damping, seed, t_sim_relax, args.grid_n,
                            pbc=bool(args.pbc), use_gpu=bool(args.use_gpu),
                            damping_noise=args.damping_noise,
                            needle_frac=args.needle_frac, needle_thickness=args.needle_thickness,
                            steel_density=args.steel_density, steel_Bsat=args.steel_Bsat,
                            moment=args.moment, inertia=args.inertia,
                            dt_factor=args.dt_factor, noise=args.noise,
                            domain_tol=args.domain_tol,
                            pivot_radius=args.pivot_radius, pivot_thickness=args.pivot_thickness,
                            pivot_density=args.pivot_density, pivot_mass=args.pivot_mass)

                    csv_path = os.path.join(relax_dir, tag + '.csv')
                    meta_path = os.path.join(meta_dir, tag + '_relax.json')
                    save_field_log_csv(field_log, csv_path)
                    with open(meta_path, 'w') as f:
                        json.dump(meta, f, indent=2)

                    manifest_rows.append(dict(
                        tag=tag, stage='relax', geometry=geometry,
                        damp_idx=damp_idx, seed=seed, damping=damping,
                        Q=meta['Q'], n_steps=meta['n_steps'],
                        wall_seconds=meta['wall_seconds'],
                        csv_path=csv_path, meta_path=meta_path,
                    ))
                    _bar(run_idx2, total_runs_relax, label="RELAX",
                         extra=f"{tag}  Q={meta['Q']:.3f}  nd={meta['n_domains']}")

    # ── Stage 3: FORC measurement ──────────────────────────────────
    if not args.skip_forc_stage:
        total_runs_forc = len(geometries) * len(damping_vals) * args.n_seeds
        run_idx3 = 0
        _bar(0, total_runs_forc, label="FORC")
        for geometry in geometries:
            for damp_idx, damping in enumerate(damping_vals):
                for seed in range(args.n_seeds):
                    run_idx3 += 1
                    tag = f"{geometry}_damp{damp_idx:02d}_seed{seed:02d}"

                    with open(os.devnull, 'w') as _dn, \
                            contextlib.redirect_stdout(_dn):
                        field_log, meta = run_one_forc(
                            geometry, args.R, damping, seed, B_max, args.grid_n,
                            pbc=bool(args.pbc), use_gpu=bool(args.use_gpu),
                            damping_noise=args.damping_noise,
                            needle_frac=args.needle_frac, needle_thickness=args.needle_thickness,
                            steel_density=args.steel_density, steel_Bsat=args.steel_Bsat,
                            moment=args.moment, inertia=args.inertia,
                            dt_factor=args.dt_factor, noise=args.noise,
                            pivot_radius=args.pivot_radius, pivot_thickness=args.pivot_thickness,
                            pivot_density=args.pivot_density, pivot_mass=args.pivot_mass,
                            forc_n_curves=args.forc_n_curves,
                            forc_Br_min=args.forc_Br_min,
                            forc_rate=args.forc_rate,
                            forc_t_sat=args.forc_t_sat,
                            forc_t_ramp_down=args.forc_t_ramp_down,
                            forc_t_ramp_up=args.forc_t_ramp_up)

                    csv_path = os.path.join(forc_dir, tag + '.csv')
                    meta_path = os.path.join(meta_dir, tag + '_forc.json')
                    save_forc_csv(field_log, csv_path)
                    with open(meta_path, 'w') as f:
                        json.dump(meta, f, indent=2)

                    manifest_rows.append(dict(
                        tag=tag, stage='forc', geometry=geometry,
                        damp_idx=damp_idx, seed=seed, damping=damping,
                        Q=meta['Q'], n_steps=meta['n_steps'],
                        wall_seconds=meta['wall_seconds'],
                        csv_path=csv_path, meta_path=meta_path,
                    ))
                    _bar(run_idx3, total_runs_forc, label="FORC",
                         extra=f"{tag}  Q={meta['Q']:.3f}  {meta['wall_seconds']:.1f}s")

    # ── manifest.csv ──────────────────────────────────────────────────
    import csv as csv_module
    if manifest_rows:
        keys = list(manifest_rows[0].keys())
        with open(manifest_path, 'w', newline='') as f:
            writer = csv_module.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            for row in manifest_rows:
                writer.writerow(row)

    total_wall = time.perf_counter() - t_campaign_start
    print()
    print(f"Campaign completed in {total_wall/60.0:.1f} min.")
    print(f"Manifest: {manifest_path}")


if __name__ == '__main__':
    main()
