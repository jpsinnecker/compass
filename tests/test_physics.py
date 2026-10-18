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


def test_quality_factor_is_independent_of_field_amplitude(tmp_path):
    common = ["--geometry", "square", "--N", "4", "--M", "4", "--field_mode", "static",
              "--damping", "1e-7", "--t_sim", "0.01"]
    _, meta_lo = _run(tmp_path, "q_lo", *common, "--B_max_factor", "1")
    _, meta_hi = _run(tmp_path, "q_hi", *common, "--B_max_factor", "8")
    assert meta_lo["derived"]["Q"] == pytest.approx(meta_hi["derived"]["Q"], rel=1e-9)
    # ... and on the protocol itself
    common_forc = [a for a in common if a not in ("static",)]
    common_forc[common_forc.index("--field_mode") + 1:common_forc.index("--field_mode") + 1] = ["forc"]
    _, meta_forc = _run(tmp_path, "q_forc", *common_forc, "--B_max_factor", "8",
                        "--forc_n_curves", "2", "--forc_t_sat", "0.002",
                        "--forc_t_ramp_down", "0.002", "--forc_t_ramp_up", "0.002")
    assert meta_forc["derived"]["Q"] == pytest.approx(meta_lo["derived"]["Q"], rel=1e-9)


def test_quality_factor_uses_lattice_field_scale(tmp_path):
    # Q0 = sqrt(m B_lat I) / b with B_lat the intrinsic lattice field; for a
    # 20x20 square array B_lat ~ 2 B_ref, so Q0 / Q_ref = sqrt(B_lat / B_ref).
    _, meta = _run(tmp_path, "q_lat", "--geometry", "square", "--N", "20", "--M", "20",
                   "--field_mode", "static", "--damping", "1e-7", "--t_sim", "0.005")
    d = meta["derived"]
    ratio_B = d["B_lat_T"] / d["B_ref_T"]
    assert 1.8 < ratio_B < 2.3
    assert d["Q"] / d["Q_ref"] == pytest.approx(math.sqrt(ratio_B), rel=1e-9)


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


def _balance_residual(tmp_path, tag, dt_factor):
    _, meta = _run(
        tmp_path, tag,
        "--geometry", "triangular", "--N", "5", "--M", "5",
        "--field_mode", "sine", "--B_max_factor", "4", "--field_freq", "3",
        "--damping", "5e-8", "--t_sim", "0.5", "--dt_factor", str(dt_factor),
    )
    return meta["energy_ledger"]["balance_residual_rel"]


def test_energy_balance_residual_is_second_order(tmp_path):
    coarse = _balance_residual(tmp_path, "bal_coarse", 0.02)
    fine = _balance_residual(tmp_path, "bal_fine", 0.01)
    assert fine < 1e-3
    assert 2.5 < coarse / fine < 6.0, f"residual ratio {coarse / fine:.2f}, expected ~4"


def test_hysteresis_loop_area_equals_dissipated_energy(tmp_path):
    # Over a closed field cycle from a settled state the field work equals the
    # dissipation (up to the change of internal energy), i.e. the loop area
    # -K m * closed integral of M_proj dB equals E_diss over that cycle.
    rows, meta = _run(
        tmp_path, "loop",
        "--geometry", "square", "--N", "5", "--M", "5",
        "--field_mode", "hysteresis", "--B_max_factor", "4", "--damping", "2e-7",
        "--t_sim", "2.0", "--log_every", "1", "--dt_factor", "0.01",
    )
    branch = [r["branch"] for r in rows]
    # closed loop: start of "pos_to_0" (at +Bmax) to the end of "0_to_pos_final" (+Bmax again)
    i0 = branch.index("pos_to_0") - 1
    i1 = len(rows) - 1
    B = _col(rows, "B_scalar_T")[i0:i1 + 1]
    M = _col(rows, "M_proj")[i0:i1 + 1]
    W = _col(rows, "W_field_J")
    D = _col(rows, "E_diss_J")
    E = _col(rows, "E_total_J")
    K = meta["derived"]["K"]
    moment = meta["config"]["moment"]
    area_work = -K * moment * np.trapz(M, B)
    assert area_work == pytest.approx(W[i1] - W[i0], rel=2e-2)
    assert (W[i1] - W[i0]) - (D[i1] - D[i0]) == pytest.approx(E[i1] - E[i0], abs=1e-3 * (D[i1] - D[i0]))



# ----------------------------------------------------------------------------
# dimensionless sweep rate r = (dB/dt) tau_relax / B_ref (roadmap item 4)
# ----------------------------------------------------------------------------

def _closed_loop_area_and_work(rows, meta):
    """Loop area -K m * closed integral of M_proj dB and field work over the
    closed cycle +Bmax -> -Bmax -> +Bmax (start of pos_to_0 to the end)."""
    branch = [r["branch"] for r in rows]
    i0 = branch.index("pos_to_0") - 1
    i1 = len(rows) - 1
    B = _col(rows, "B_scalar_T")[i0:i1 + 1]
    M = _col(rows, "M_proj")[i0:i1 + 1]
    W = _col(rows, "W_field_J")
    area = -meta["derived"]["K"] * meta["config"]["moment"] * np.trapz(M, B)
    return area, W[i1] - W[i0]


def test_sweep_rate_dimless_sets_ramp_and_is_recorded(tmp_path):
    common = ["--geometry", "square", "--N", "4", "--M", "4", "--field_mode", "hysteresis",
              "--B_max_factor", "3", "--damping", "1e-7", "--log_every", "50"]
    _, meta = _run(tmp_path, "r_in", *common, "--sweep_rate_dimless", "0.5")
    sr, d = meta["sweep_rate"], meta["derived"]
    assert sr["input"] == "sweep_rate_dimless"
    assert sr["sweep_rate_dimless"] == pytest.approx(0.5, rel=1e-12)
    tau = max(1.0 / d["omega_lat_rad_s"], 1e-7 / (meta["config"]["moment"] * d["B_lat_T"]))
    assert sr["tau_relax_s"] == pytest.approx(tau, rel=1e-12)
    B_ext = meta["config"]["B_ext"]
    assert meta["config"]["t_sim"] == pytest.approx(5.0 * B_ext * tau / (0.5 * d["B_ref_T"]), rel=1e-12)

    # Recorded also when the ramp comes from --t_sim, and round-trips.
    _, meta_t = _run(tmp_path, "r_t", *common, "--t_sim", str(meta["config"]["t_sim"]))
    assert meta_t["sweep_rate"]["input"] == "t_sim"
    assert meta_t["sweep_rate"]["sweep_rate_dimless"] == pytest.approx(0.5, rel=1e-9)

    # ... and for FORC, with the flag and with fixed ramp times.
    forc = ["--geometry", "square", "--N", "4", "--M", "4", "--field_mode", "forc",
            "--B_max_factor", "3", "--damping", "1e-7", "--forc_n_curves", "2",
            "--forc_t_sat", "0.002", "--log_every", "50"]
    _, meta_f = _run(tmp_path, "r_forc", *forc, "--sweep_rate_dimless", "2")
    assert meta_f["sweep_rate"]["sweep_rate_dimless"] == pytest.approx(2.0, rel=1e-12)
    assert meta_f["derived"]["sweep_rate_T_per_s"] == pytest.approx(meta_f["sweep_rate"]["sweep_rate_T_per_s"])
    _, meta_ff = _run(tmp_path, "r_forc_fixed", *forc, "--forc_t_ramp_down", "0.01", "--forc_t_ramp_up", "0.01")
    assert meta_ff["sweep_rate"]["input"] == "forc_ramp_times"
    assert meta_ff["sweep_rate"]["sweep_rate_dimless"] > 0.0

    with pytest.raises(SystemExit):
        _run(tmp_path, "r_both", *common, "--sweep_rate_dimless", "0.5", "--t_sim", "1")
    with pytest.raises(SystemExit):
        _run(tmp_path, "r_forc_both", *forc, "--sweep_rate_dimless", "2", "--forc_rate", "0.1")


def test_loop_area_converges_in_the_quasi_static_limit(tmp_path):
    # At fixed Q0 the loop area must approach a finite quasi-static value
    # A0 = lim_{r->0} A(r): successive differences over decades of r shrink.
    # 6x6 triangular, Q0 = 2. Measured (Oct 2026): A = 1.03e-3, 2.80e-4,
    # 1.62e-4 J at r = 0.1, 0.01, 0.001 (and 1.46e-4 at r = 3e-4), dt-converged
    # to 1e-5 at dt_factor 0.05. The closed-cycle field work must equal the area.
    geo = ["--geometry", "triangular", "--N", "6", "--M", "6"]
    _, probe = _run(tmp_path, "qs_probe", *geo, "--field_mode", "static", "--t_sim", "0.001")
    Q0 = 2.0
    damping = probe["derived"]["omega_lat_rad_s"] * probe["config"]["inertia"] / Q0

    areas = []
    for r in (0.1, 0.01, 0.001):
        rows, meta = _run(
            tmp_path, f"qs_{r:g}", *geo,
            "--field_mode", "hysteresis", "--B_max_factor", "3", "--damping", repr(damping),
            "--sweep_rate_dimless", repr(r), "--dt_factor", "0.05", "--log_every", "4",
        )
        assert meta["derived"]["Q"] == pytest.approx(Q0, rel=1e-9)
        area, work = _closed_loop_area_and_work(rows, meta)
        assert work == pytest.approx(area, rel=1e-3)
        areas.append(area)

    a_fast, a_mid, a_slow = areas
    assert a_fast > a_mid > a_slow > 0.0
    assert abs(a_mid - a_slow) < 0.3 * abs(a_fast - a_mid), f"no convergence: areas {areas}"


# ----------------------------------------------------------------------------
# event-driven quasi-static loop (--field_mode hysteresis_qs)
# ----------------------------------------------------------------------------

_QS_N = 3
_QS_COMMON = ["--geometry", "square", "--N", str(_QS_N), "--M", str(_QS_N), "--field_mode", "hysteresis_qs",
              "--B_max_factor", "2.5", "--qs_dB_factor", "0.05", "--qs_tol_torque", "1e-4", "--qs_tol_omega", "1e-4"]


def _damping_for_Q0(Q0, N=_QS_N, geometry="square"):
    geom = _lattice(N, geometry)
    backend = compass.select_backend(False)
    tensor = _tensor(geom, FULL_CUTOFF_SHELLS)
    moment = compass.compute_moment_from_geometry(geom.needle_len, geom.needle_width, compass.NEEDLE_THICKNESS_DEFAULT)
    inertia = compass.compute_inertia_from_geometry(geom.needle_len, geom.needle_width, compass.NEEDLE_THICKNESS_DEFAULT)
    B_lat = moment * compass.lattice_field_per_moment(tensor, backend)
    return math.sqrt(moment * B_lat * inertia) / Q0


def _qs_area(rows, meta):
    i_top = max(i for i, r in enumerate(rows) if r["branch"] == "0_to_pos")
    W = _col(rows, "W_field_J")
    norm = meta["derived"]["K"] * meta["config"]["moment"] * meta["derived"]["B_ref_T"]
    return (W[-1] - W[i_top]) / norm


def test_qs_path_shape():
    path = compass.build_qs_path(3.0, 1.0)
    assert [b for b, _ in path] == [1, 2, 3, 2, 1, 0, -1, -2, -3, -2, -1, 0, 1, 2, 3]
    assert path[0][1] == "0_to_pos" and path[-1][1] == "0_to_pos_final"


def test_qs_logs_rest_states_that_are_equilibria(tmp_path):
    rows, meta = _run(tmp_path, "qs_eq", *_QS_COMMON, "--damping", repr(_damping_for_Q0(3.0)))
    qs = meta["quasi_static"]
    assert qs["unconverged_field_steps"] == 0
    assert len(rows) == qs["n_field_steps"] + 1            # initial row + one per field step
    # The final state, re-evaluated independently, is a mechanical equilibrium.
    st = np.load(tmp_path / "states" / "qs_eq_final.npz")
    geom = _lattice(_QS_N)
    Bx, By = _field(_tensor(geom, FULL_CUTOFF_SHELLS), meta["config"]["moment"], st["theta"])
    bx = float(rows[-1]["Bx_T"]); by = float(rows[-1]["By_T"])
    m = meta["config"]["moment"]
    torque = m * np.cos(st["theta"]) * (By + by) - m * np.sin(st["theta"]) * (Bx + bx)
    assert np.max(np.abs(torque)) <= qs["tol_torque_Nm"] * 1.0001
    assert np.max(np.abs(st["omega"])) <= qs["tol_omega_rad_s"] * 1.0001


def test_qs_energy_balance_includes_field_jumps(tmp_path):
    _, meta = _run(tmp_path, "qs_bal", *_QS_COMMON, "--damping", repr(_damping_for_Q0(3.0)))
    assert meta["energy_ledger"]["balance_residual_rel"] < 2e-3


@pytest.mark.slow
def test_qs_overdamped_limit_matches_low_Q0_inertial(tmp_path):
    # In the quasi-static limit, inertial dynamics at small Q0 and the gradient
    # flow (Q0 -> 0) select the same sequence of minima.
    # Needs the validated rest tolerance: a loose one stops early near
    # instabilities, and differently for the two rest criteria.
    b = repr(_damping_for_Q0(0.2))
    tight = ["--qs_tol_torque", "1e-5", "--qs_tol_omega", "1e-5"]
    rows_i, meta_i = _run(tmp_path, "qs_in", *_QS_COMMON, *tight, "--damping", b)
    rows_o, meta_o = _run(tmp_path, "qs_od", *_QS_COMMON, *tight, "--damping", b, "--qs_dynamics", "overdamped")
    assert meta_o["quasi_static"]["unconverged_field_steps"] == 0
    a_i, a_o = _qs_area(rows_i, meta_i), _qs_area(rows_o, meta_o)
    assert a_o == pytest.approx(a_i, rel=0.03)
    M_i, M_o = _col(rows_i, "M_proj"), _col(rows_o, "M_proj")
    assert np.max(np.abs(M_i - M_o)) < 0.05


def test_qs_rejects_sweep_rate(tmp_path):
    with pytest.raises(SystemExit):
        _run(tmp_path, "qs_bad", *_QS_COMMON, "--sweep_rate_dimless", "0.01")


def test_qs_state_log_matches_csv_and_stop_after(tmp_path):
    # --qs_state_log stores one rest state per CSV row; --qs_stop_after
    # truncates the loop without changing the rows that are computed.
    b = repr(_damping_for_Q0(3.0))
    rows_full, _ = _run(tmp_path, "qs_full", *_QS_COMMON, "--damping", b)
    rows, meta = _run(tmp_path, "qs_log", *_QS_COMMON, "--damping", b, "--qs_state_log", "--qs_stop_after", "7")
    assert len(rows) == 8 and meta["quasi_static"]["stopped_early"]
    assert [r["M_proj"] for r in rows] == [r["M_proj"] for r in rows_full[:8]]
    st = np.load(tmp_path / "states" / "qs_log_qs_states.npz")
    assert st["theta"].shape == (8, meta["derived"]["K"])
    M = np.cos(st["theta"] - float(st["phi_ext"])).mean(axis=1)
    np.testing.assert_allclose(M, _col(rows, "M_proj"), atol=1e-10)
    np.testing.assert_array_equal(st["theta"][-1], np.load(tmp_path / "states" / "qs_log_final.npz")["theta"])
