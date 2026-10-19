"""Bench analysis (bench_fit.py): each estimator recovers known parameters from
synthetic data generated with the same equation of motion."""
import math

import numpy as np
import pytest

import bench_fit as bf

I, M, B = 1.63e-10, 4.9e-3, 1e-4          # V79 needle in a 100 uT coil field
MB = M * B


def test_geometry_and_added_mass():
    # Rhombus plate: I = M (L^2 + W^2) / 24.
    assert bf.inertia_rhombus(0.1e-3, 0.01, 0.003) == pytest.approx(0.1e-3 * (1e-4 + 9e-6) / 24)
    # Added inertia: omega^2 = mB/I, so I/(I + dI) = (f1/f0)^2.
    dI = 0.4 * I
    f0 = math.sqrt(MB / I) / (2 * math.pi)
    f1 = math.sqrt(MB / (I + dI)) / (2 * math.pi)
    assert bf.inertia_added_mass(f0, f1, dI) == pytest.approx(I, rel=1e-12)


def test_freq_vs_field_recovers_m_over_I_and_residual_field():
    B_coil = np.linspace(2e-5, 3e-4, 7)
    B_res = 1.8e-5                          # uncompensated horizontal Earth field
    f = np.sqrt(M / I * (B_coil + B_res)) / (2 * math.pi)
    k, b_res, rel = bf.fit_freq_vs_field(B_coil, f)
    assert k == pytest.approx(M / I, rel=1e-10)
    assert b_res == pytest.approx(B_res, rel=1e-8)
    assert rel < 1e-12


def test_dipole_field_recovers_moment():
    r = np.linspace(0.03, 0.08, 6)
    Bz = 2e-7 * M / r**3 + 3e-6
    m, bg = bf.fit_dipole(r, Bz)
    assert m == pytest.approx(M, rel=1e-10)
    assert bg == pytest.approx(3e-6, rel=1e-6)


def test_deadband_recovers_static_friction():
    tau_s = 8e-9
    half = math.asin(tau_s / MB)
    rng = np.random.default_rng(1)
    th = 0.3 + rng.uniform(-half, half, 40)
    est, th0, _ = bf.tau_static_from_deadband(th, M, B)
    assert est == pytest.approx(tau_s, rel=0.1)
    assert th0 == pytest.approx(0.3, abs=0.1 * half)


@pytest.mark.parametrize("noise_deg", [0.0, 0.2])
def test_decay_fit_recovers_viscous_and_coulomb_terms(noise_deg):
    b_true, tc_true = 2e-10, 5e-9           # Q of one needle ~45; tau_c/(mB) ~ 1 %
    t, th = bf.simulate(I, M, B, b_true, tc_true, 8e-9, math.radians(30), 3.0, 240.0,
                        math.radians(noise_deg), seed=3)
    w0 = math.sqrt(MB / I)
    ths = bf.smooth(t, th, 2 * math.pi / w0)
    noise = float(np.std(th - ths))
    n = bf.stuck_index(t, ths, w0, noise=noise)
    tk, ak, _ = bf.find_extrema(t[:n], ths[:n], min_sep=0.25 * 2 * math.pi / w0)
    b_env, tc_env, T = bf.fit_envelope(tk, ak, I, MB, a_min=5 * noise)
    # The envelope (small-amplitude map) is only a starting point: within 25 %.
    assert b_env == pytest.approx(b_true, rel=0.25)
    assert tc_env == pytest.approx(tc_true, rel=0.25)
    fit = bf.fit_decay(t[:n], th[:n], I, MB, b_env, tc_env)
    assert fit["b"] == pytest.approx(b_true, rel=0.04)
    assert fit["tau_c"] == pytest.approx(tc_true, rel=0.04)
