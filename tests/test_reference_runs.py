"""Regression tests: re-run each small config.yaml-variant case and compare
against the checked-in reference data in tests/reference/.

These are NOT physics-correctness tests -- the cases use large timesteps
deliberately, purely for speed (see tests/cases.py and tests/README.md), and
in fact trip compassV2_2.py's own numerical-stability monitor (a printed
warning, not a failure). What's being tested is that the same config.yaml
variant + the same code always produces the same numbers: a config or
engine change that alters simulation output will fail these tests, which is
the point.

If a failure here is an *intentional* behavior change, regenerate the
references (see tests/README.md) rather than hand-editing them.
"""

import pytest

from cases import CASES
from helpers import FIXTURES_DIR, REFERENCE_DIR, compare_csv, compare_metadata_derived, compare_npz_state, run_case


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_case_matches_reference(case, tmp_path):
    name = case["name"]
    fixture_path = FIXTURES_DIR / f"config_{name}.yaml"
    ref_dir = REFERENCE_DIR / name

    assert fixture_path.exists(), (
        f"missing fixture {fixture_path} -- run `python3 tests/generate_reference.py` first"
    )
    assert ref_dir.exists(), (
        f"missing reference dir {ref_dir} -- run `python3 tests/generate_reference.py` first"
    )

    csv_path, meta_path, state_path = run_case(fixture_path, tmp_path)

    expected_csv = ref_dir / "data" / csv_path.name
    expected_meta = ref_dir / "meta" / meta_path.name
    expected_state = ref_dir / "states" / state_path.name

    for expected in (expected_csv, expected_meta, expected_state):
        assert expected.exists(), f"reference file missing: {expected}"

    compare_csv(csv_path, expected_csv)
    compare_metadata_derived(meta_path, expected_meta)
    compare_npz_state(state_path, expected_state)
