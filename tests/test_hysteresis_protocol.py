"""tests/test_hysteresis_protocol.py

Fast, protocol-level unit tests for the hysteresis field schedule
(build_hysteresis_schedule / FieldProtocol.at, mode="hysteresis").

These do NOT run a full simulation: they call the schedule builder and the
protocol's .at(t) sampling function directly, which is enough to verify the
three properties the feature is required to have:

  1. With no slow window (the default), the new schedule-based FieldProtocol
     reproduces the previous fixed-five-equal-segment formula exactly
     (bit-for-bit), so existing datasets and behavior are unaffected.
  2. With a slow window, the sweep rate strictly inside the window equals
     base_rate / hyst_slow_factor exactly.
  3. B(t) is continuous everywhere, including across window boundaries
     (dB/dt is allowed to jump there, but B itself must not).

See FORC.md's "variable sweep rate trap" (V67-V69): sweep rates must be
explicit, constant within each declared segment, and reported in metadata --
never an implicit byproduct of holding total time fixed. The tests below
check exactly that the reported per-segment rate matches the rate actually
used to advance B(t).
"""

from __future__ import annotations

import math
from pathlib import Path

import helpers

REPO_ROOT = Path(__file__).resolve().parent.parent


def _legacy_hysteresis_B(t: float, Bmax: float, t_sim: float) -> float:
    """Independent re-implementation of the PRE-EXISTING fixed-five-equal-
    segment hysteresis formula (linear spacing), kept deliberately separate
    from compassV2_2.py's source so this test does not just check the code
    against itself."""
    T = max(t_sim, 1e-30)
    t5 = T / 5.0
    if t <= t5:
        u, sgn = t / t5, +1.0
    elif t <= 2.0 * t5:
        u, sgn = 1.0 - (t - t5) / t5, +1.0
    elif t <= 3.0 * t5:
        u, sgn = (t - 2.0 * t5) / t5, -1.0
    elif t <= 4.0 * t5:
        u, sgn = 1.0 - (t - 3.0 * t5) / t5, -1.0
    else:
        u, sgn = (t - 4.0 * t5) / t5, +1.0
    u = max(0.0, min(1.0, u))
    return sgn * Bmax * u


def _load_module():
    return helpers.load_fresh_compassV2_2(helpers.BASE_CONFIG_PATH)


def test_no_window_bit_identical_to_legacy_formula():
    mod = _load_module()
    Bmax = 0.0123
    t_sim = 0.77

    schedule = mod.build_hysteresis_schedule(Bmax, t_sim, None, None, 1.0)
    protocol = mod.FieldProtocol(
        mode="hysteresis", Bmax=Bmax, phi=0.0, t_sim=schedule.total_time,
        hyst_schedule=schedule,
    )

    assert schedule.total_time == t_sim  # exact, not just close

    n = 4001
    for i in range(n):
        t = t_sim * i / (n - 1)
        expected = _legacy_hysteresis_B(t, Bmax, t_sim)
        actual = protocol.at(t).B_scalar
        assert actual == expected, f"mismatch at t={t}: {actual!r} != {expected!r}"

    # And past the nominal end of the cycle (old code's final "else" branch
    # clamps u to 1.0, saturating at +Bmax).
    for t in (t_sim, t_sim * 1.3, t_sim * 5.0):
        expected = _legacy_hysteresis_B(t, Bmax, t_sim)
        actual = protocol.at(t).B_scalar
        assert actual == expected, f"mismatch at t={t}: {actual!r} != {expected!r}"


def test_slow_window_rate_equals_base_rate_over_factor():
    mod = _load_module()
    Bmax = 0.05
    t_sim_nominal = 1.0
    factor = 6.0
    lo, hi = 0.01, 0.03

    schedule = mod.build_hysteresis_schedule(Bmax, t_sim_nominal, lo, hi, factor)
    base_rate = schedule.base_rate_T_per_s
    assert base_rate == Bmax / (t_sim_nominal / 5.0)
    assert schedule.slow_rate_T_per_s == base_rate / factor

    # leg.rate_T_per_s is recomputed from the leg's own (cumulative-sum)
    # t_start/t_end, so it can differ from base_rate/factor by a few ULPs
    # even though both express the same intended rate; compare with a tight
    # relative tolerance rather than bit-exact equality (consistent with
    # tests/helpers.py's RTOL convention for this repo's float comparisons).
    found_inside = False
    for leg in schedule.legs:
        mag_mid = 0.5 * abs(leg.B_start + leg.B_end)
        if lo < mag_mid < hi:
            found_inside = True
            assert math.isclose(abs(leg.rate_T_per_s), base_rate / factor, rel_tol=1e-9)
        else:
            assert math.isclose(abs(leg.rate_T_per_s), base_rate, rel_tol=1e-9)
    assert found_inside, "no leg fell inside the slow window; test setup is wrong"

    # total simulated time must grow to accommodate the slowed segments.
    assert schedule.total_time > t_sim_nominal


def test_B_continuous_across_window_boundaries():
    mod = _load_module()
    Bmax = 0.02
    t_sim_nominal = 0.5
    factor = 10.0
    lo, hi = 0.004, 0.012

    schedule = mod.build_hysteresis_schedule(Bmax, t_sim_nominal, lo, hi, factor)
    protocol = mod.FieldProtocol(
        mode="hysteresis", Bmax=Bmax, phi=0.0, t_sim=schedule.total_time,
        hyst_schedule=schedule,
    )

    eps = 1e-9
    for leg in schedule.legs:
        b_end = protocol.at(leg.t_end).B_scalar
        b_just_before = protocol.at(leg.t_end - eps).B_scalar
        b_just_after = protocol.at(leg.t_end + eps).B_scalar
        assert math.isclose(b_end, b_just_before, abs_tol=1e-6)
        assert math.isclose(b_end, b_just_after, abs_tol=1e-6)

    # Dense scan across the whole cycle: no sample-to-sample jump should
    # exceed what a continuous piecewise-linear B(t) can produce over a
    # small dt (i.e., no discontinuity anywhere).
    n = 5000
    T = schedule.total_time
    prev = protocol.at(0.0).B_scalar
    max_step = 2.0 * Bmax * (T / n) / max(schedule.slow_rate_T_per_s and (T / 5.0) or 1.0, 1.0)
    for i in range(1, n + 1):
        t = T * i / n
        cur = protocol.at(t).B_scalar
        assert abs(cur - prev) < 0.05 * Bmax  # generous bound; catches real jumps
        prev = cur


def test_hyst_slow_window_parsing():
    mod = _load_module()
    assert mod._parse_hyst_slow_window(None) == (None, None)
    assert mod._parse_hyst_slow_window("0.01,0.03") == (0.01, 0.03)
    try:
        mod._parse_hyst_slow_window("bad")
        assert False, "expected SystemExit for malformed spec"
    except SystemExit:
        pass
    try:
        mod._parse_hyst_slow_window("0.03,0.01")
        assert False, "expected SystemExit for B_lo >= B_hi"
    except SystemExit:
        pass
