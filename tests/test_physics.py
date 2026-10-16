"""tests/test_physics.py

Physical-validity tests for compass.py (unlike test_reference_runs.py, which
only checks that the numbers did not change).

Two groups:

* Passing tests: invariants the engine must keep satisfying (tensor == direct
  pair sum, torque == -dE/dtheta, symplectic energy conservation at b=0 with
  second-order convergence, monotonic energy decay with damping).

* xfail(strict=True) tests: known physics defects documented in
  docs/AUDIT_PHYSICS_ROADMAP.md (sections 2.2-2.4 and 2.6). Each one encodes
  the *target* behaviour. When a Phase-1 fix lands, the corresponding test
  starts passing, strict xfail turns that into an XPASS failure, and the
  marker must be removed in the same commit. That keeps the list of known
  defects honest.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
import pytest

import helpers

compass = helpers.load_fresh_compass(helpers.BASE_CONFIG_PATH)

MU0_OVER_4PI = 1.0e-7
R = 0.0065                    # default half pivot spacing [m]
FULL_CUTOFF_SHELLS = 1.0e3    # explicit, effectively infinite cutoff (independent of the CLI default)


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------

def _lattice(n, geometry="square", noise=0.0, seed=0):
    rng = np.random.default_rng(seed)
    return compass.make_lattice(n, n, geometry, R, 0.8, noise, rng)


def _tensor(geom, cutoff_shells):
    backend = compass.select_backend(False)
    return compass.precompute_dipolar_tensor(
        geom.xs, geom.ys, cutoff_shells * geom.r_nn, cutoff_shells, backend
    )


def _field(tensor, moment, theta):
    mx = moment * np.cos(theta)
    my = moment * np.sin(theta)
    return tensor.Axx @ mx + tensor.Axy @ my, tensor.Axy @ mx + tensor.Ayy @ my


def _run(tmp_path, tag, *extra):
    """Run compass.main() in tmp_path and return (rows, metadata)."""
    argv = [
        "--out_dir", str(tmp_path), "--tag", tag,
        "--png_dpi", "20", "--seed", "11", "--damping_noise", "0",
        *extra,
    ]
    compass.main(argv)
    with open(tmp_path / "data" / f"{tag}.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    with open(tmp_path / "meta" / f"{tag}.json") as fh:
        meta = json.load(fh)
    return rows, meta


def _col(rows, name):
    return np.array([float(r[name]) for r in rows])


def _max_energy_drift(tmp_path, tag, dt_factor):
    rows, _ = _run(
        tmp_path, tag,
        "--geometry", "square", "--N", "6", "--M", "6",
        "--field_mode", "static", "--B_ext", "0", "--damping", "0",
        "--noise", "1.0", "--t_sim", "0.3", "--log_every", "5",
        "--dt_factor", str(dt_factor), "--cutoff_shells", str(FULL_CUTOFF_SHELLS),
    )
    E = _col(rows, "E_total_J")
    scale = abs(_col(rows, "E_dip_J")[0])
    return np.max(np.abs(E - E[0])) / scale


# ----------------------------------------------------------------------------
# invariants that must hold
# ----------------------------------------------------------------------------

def test_tensor_matches_direct_pair_sum():
    geom = _lattice(5, "triangular")
    moment = 1.0e-3
    theta = np.random.default_rng(1).uniform(-math.pi, math.pi, geom.K)
    Bx, By = _field(_tensor(geom, FULL_CUTOFF_SHELLS), moment, theta)

    m = moment * np.stack([np.cos(theta), np.sin(theta)], axis=1)
    pos = np.stack([geom.xs, geom.ys], axis=1)
    B_direct = np.zeros_like(m)
    for i in range(geom.K):
        for j in range(geom.K):
            if i == j:
                continue
            r = pos[i] - pos[j]
            d = np.linalg.norm(r)
            rh = r / d
            B_direct[i] += MU0_OVER_4PI * (3.0 * np.dot(m[j], rh) * rh - m[j]) / d**3

    np.testing.assert_allclose(Bx, B_direct[:, 0], rtol=1e-12, atol=1e-18)
    np.testing.assert_allclose(By, B_direct[:, 1], rtol=1e-12, atol=1e-18)


def test_torque_is_minus_energy_gradient():
    geom = _lattice(4, "honeycomb")
    moment = 1.0e-3
    bx, by = 2.0e-4, -1.0e-4
    tensor = _tensor(geom, FULL_CUTOFF_SHELLS)
    theta = np.random.default_rng(2).uniform(-math.pi, math.pi, geom.K)

    def energy(th):
        Bx, By = _field(tensor, moment, th)
        mx, my = moment * np.cos(th), moment * np.sin(th)
        return -0.5 * np.sum(mx * Bx + my * By) - np.sum(mx * bx + my * by)

    Bx, By = _field(tensor, moment, theta)
    torque = moment * np.cos(theta) * (By + by) - moment * np.sin(theta) * (Bx + bx)

    h = 1e-6
    grad = np.empty(geom.K)
    for k in range(geom.K):
        tp, tm = theta.copy(), theta.copy()
        tp[k] += h
        tm[k] -= h
        grad[k] = (energy(tp) - energy(tm)) / (2 * h)

    np.testing.assert_allclose(torque, -grad, rtol=1e-6, atol=1e-6 * np.max(np.abs(torque)))


def test_energy_conservation_at_zero_damping_is_second_order(tmp_path):
    drift_coarse = _max_energy_drift(tmp_path, "e_coarse", 0.01)
    drift_fine = _max_energy_drift(tmp_path, "e_fine", 0.005)
    assert drift_fine < 2e-3, f"energy drift {drift_fine:.2e} too large for dt_factor=0.005"
    ratio = drift_coarse / drift_fine
    assert 2.5 < ratio < 6.0, f"drift ratio {ratio:.2f}, expected ~4 for a 2nd-order symplectic scheme"


def test_damped_relaxation_dissipates_energy(tmp_path):
    rows, _ = _run(
        tmp_path, "damped",
        "--geometry", "square", "--N", "6", "--M", "6",
        "--field_mode", "static", "--B_ext", "0", "--damping", "2e-7",
        "--noise", "1.0", "--t_sim", "0.5", "--log_every", "5",
        "--dt_factor", "0.005", "--cutoff_shells", str(FULL_CUTOFF_SHELLS),
    )
    E = _col(rows, "E_total_J")
    scale = abs(_col(rows, "E_dip_J")[0])
    assert E[-1] < E[0]
    assert np.max(np.diff(E)) < 1e-4 * scale, "total energy increased during free damped relaxation"


# ----------------------------------------------------------------------------
# known defects (docs/AUDIT_PHYSICS_ROADMAP.md) -- remove the marker when fixed
# ----------------------------------------------------------------------------

def test_default_cutoff_reproduces_full_dipolar_sum():
    geom = _lattice(20)
    moment = 1.0e-3
    default_shells = compass.build_parser().parse_args([]).cutoff_shells
    if default_shells is None or default_shells <= 0.0:   # documented sentinel: full pair sum
        default_shells = math.inf
    theta = np.zeros(geom.K)  # uniform state along +x
    Bx_def, _ = _field(_tensor(geom, default_shells), moment, theta)
    Bx_full, _ = _field(_tensor(geom, FULL_CUTOFF_SHELLS), moment, theta)
    centre = np.argmin((geom.xs - geom.xs.mean()) ** 2 + (geom.ys - geom.ys.mean()) ** 2)
    assert abs(Bx_def[centre] / Bx_full[centre] - 1.0) < 0.02


def test_full_sum_default_runs_and_writes_valid_json(tmp_path):
    # Default CLI run (no cutoff flags) must use the full sum and keep the
    # metadata strict-JSON (no Infinity) so other tools can parse it.
    _, meta = _run(tmp_path, "fullsum", "--geometry", "honeycomb", "--N", "4", "--M", "4",
                   "--field_mode", "static", "--t_sim", "0.01")
    assert meta["derived"]["cutoff_policy"] == "full_sum"
    assert meta["config"]["cutoff_m"] is None
    text = (tmp_path / "meta" / "fullsum.json").read_text()
    assert "Infinity" not in text and "NaN" not in text


def test_default_timestep_resolves_local_frequency(tmp_path):
    _, meta = _run(
        tmp_path, "dt_default",
        "--geometry", "square", "--N", "6", "--M", "6",
        "--field_mode", "static", "--B_ext", "0", "--t_sim", "0.05",
    )
    assert meta["stability_monitor"]["max_omega_local_dt"] <= 0.2


def test_local_field_bound_dominates_actual_fields():
    # The dt bound must hold for any configuration: check random states and
    # the fully aligned state, on all three geometries.
    backend = compass.select_backend(False)
    moment = 1.0e-3
    rng = np.random.default_rng(5)
    for geometry in ("square", "triangular", "honeycomb"):
        geom = _lattice(6, geometry)
        tensor = _tensor(geom, FULL_CUTOFF_SHELLS)
        bound = moment * compass.local_field_bound_per_moment(tensor, backend)
        states = [rng.uniform(-math.pi, math.pi, geom.K) for _ in range(50)]
        states += [np.full(geom.K, a) for a in np.linspace(0.0, math.pi, 7)]
        for theta in states:
            Bx, By = _field(tensor, moment, theta)
            assert np.max(np.hypot(Bx, By)) <= bound * (1.0 + 1e-12)


def test_default_timestep_resolves_local_frequency_during_hysteresis(tmp_path):
    # Field sweeps reach |B_ext| = 8 B_ref; the guard must still never flag.
    _, meta = _run(
        tmp_path, "dt_hyst",
        "--geometry", "triangular", "--N", "6", "--M", "6",
        "--field_mode", "hysteresis", "--B_max_factor", "8", "--damping", "1e-8",
        "--t_sim", "0.3",
    )
    monitor = meta["stability_monitor"]
    assert monitor["flagged_steps"] == 0
    assert monitor["max_omega_local_dt"] <= 2.0 * math.pi * 0.025 + 1e-12


@pytest.mark.xfail(strict=True, reason="AUDIT_PHYSICS 2.4: Q uses B_eff=max(B_ext,B_ref) and so depends on the protocol amplitude")
def test_quality_factor_is_independent_of_field_amplitude(tmp_path):
    common = ["--geometry", "square", "--N", "4", "--M", "4", "--field_mode", "static",
              "--damping", "1e-7", "--t_sim", "0.01"]
    _, meta_lo = _run(tmp_path, "q_lo", *common, "--B_max_factor", "1")
    _, meta_hi = _run(tmp_path, "q_hi", *common, "--B_max_factor", "8")
    assert meta_lo["derived"]["Q"] == pytest.approx(meta_hi["derived"]["Q"], rel=1e-9)


@pytest.mark.xfail(strict=True, reason="AUDIT_PHYSICS 2.6: no work/dissipation ledger (W_field_J, E_diss_J) in the CSV")
def test_energy_balance_with_field_work_and_dissipation(tmp_path):
    rows, _ = _run(
        tmp_path, "balance",
        "--geometry", "square", "--N", "5", "--M", "5",
        "--field_mode", "hysteresis", "--B_max_factor", "4", "--damping", "2e-7",
        "--t_sim", "1.0", "--log_every", "2", "--dt_factor", "0.005",
        "--cutoff_shells", str(FULL_CUTOFF_SHELLS),
    )
    E = _col(rows, "E_total_J")
    W = _col(rows, "W_field_J")    # cumulative work done by the time-dependent applied field
    D = _col(rows, "E_diss_J")     # cumulative viscous dissipation, integral of b*omega^2 dt
    residual = (E - E[0]) - (W - D)
    assert np.max(np.abs(residual)) < 1e-3 * np.max(np.abs(D))
