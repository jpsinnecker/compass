"""Shared machinery for the config.yaml-variant reference tests.

Used by both generate_reference.py (offline generation) and
test_reference_runs.py (pytest re-run + compare), so the two paths cannot
silently diverge.
"""

from __future__ import annotations

import copy
import csv
import importlib.util
import json
import os
import sys
from pathlib import Path
from types import ModuleType
from typing import Any, Dict, Tuple

import numpy as np
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
BASE_CONFIG_PATH = REPO_ROOT / "config.yaml"
FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"
REFERENCE_DIR = Path(__file__).resolve().parent / "reference"

# CSV columns compared with exact equality (integers / category labels).
# Everything else in the CSV is a float compared with a numerical tolerance.
_EXACT_CSV_COLUMNS = {"step", "branch", "forc_index", "flip_field", "flip_angle"}

# Default tolerance for floating-point comparisons. The simulation is a
# deterministic fixed-seed run, so in practice values match to the bit; a
# small but nonzero tolerance is used anyway to stay robust across numpy/CPU
# combinations, as requested.
RTOL = 1e-9
ATOL = 1e-12


def deep_merge(base: Dict[str, Any], overrides: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge `overrides` into a deep copy of `base`; returns the copy."""
    result = copy.deepcopy(base)
    stack = [(result, overrides)]
    while stack:
        dst, src = stack.pop()
        for key, value in src.items():
            if isinstance(value, dict) and isinstance(dst.get(key), dict):
                stack.append((dst[key], value))
            else:
                dst[key] = copy.deepcopy(value)
    return result


def build_variant_config(overrides: Dict[str, Any], base_path: Path = BASE_CONFIG_PATH) -> Dict[str, Any]:
    """Load the repo's config.yaml and apply a case's overrides on top of it."""
    with open(base_path, "r") as fh:
        base = yaml.safe_load(fh)
    return deep_merge(base, overrides)


def write_variant_config(config_dict: Dict[str, Any], path: Path, case_name: str) -> None:
    header = (
        f"# GENERATED FIXTURE -- do not hand-edit.\n"
        f"# This is config.yaml + tests/cases.py's '{case_name}' overrides, produced by\n"
        f"# tests/generate_reference.py. Regenerate with:\n"
        f"#   python3 tests/generate_reference.py\n"
        f"# if config.yaml's schema changes or the case's overrides change.\n\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(header)
        yaml.safe_dump(config_dict, fh, sort_keys=False)


def load_fresh_compass(config_path: Path) -> ModuleType:
    """Import an independent instance of compass.py bound to `config_path`.

    compass.py (like every other script in this repo) reads config.yaml
    once at import time via `CFG = load_config()`. Setting COMPASS_CONFIG_PATH
    and loading a *fresh* module object (rather than a cached
    sys.modules['compass']) is what lets each test case run against its
    own config.yaml variant in the same process.
    """
    old_env = os.environ.get("COMPASS_CONFIG_PATH")
    os.environ["COMPASS_CONFIG_PATH"] = str(config_path)
    repo_root_str = str(REPO_ROOT)
    path_added = repo_root_str not in sys.path
    if path_added:
        sys.path.insert(0, repo_root_str)  # compass.py does `from sim_config import ...`
    try:
        mod_name = f"compass_variant_{config_path.stem}"
        spec = importlib.util.spec_from_file_location(mod_name, REPO_ROOT / "compass.py")
        module = importlib.util.module_from_spec(spec)
        # dataclasses' `from __future__ import annotations` type resolution looks
        # the module up via sys.modules[cls.__module__]; it must be registered
        # there before exec_module runs its class bodies (CPython dataclasses
        # internals, not specific to this repo's code).
        sys.modules[mod_name] = module
        try:
            spec.loader.exec_module(module)
        finally:
            del sys.modules[mod_name]
        return module
    finally:
        if path_added:
            sys.path.remove(repo_root_str)
        if old_env is None:
            os.environ.pop("COMPASS_CONFIG_PATH", None)
        else:
            os.environ["COMPASS_CONFIG_PATH"] = old_env


def run_case(config_path: Path, out_dir: Path) -> Tuple[Path, Path, Path]:
    """Run compass.py's build_parser()+run_simulation() against a config
    variant, with no CLI overrides except redirecting out_dir. Returns
    (csv_path, meta_path, state_path)."""
    module = load_fresh_compass(config_path)
    args = module.build_parser().parse_args([])
    args.out_dir = str(out_dir)
    return module.run_simulation(args)


def _read_csv_rows(path: Path):
    with open(path, "r", newline="") as fh:
        return list(csv.DictReader(fh))


def compare_csv(actual_path: Path, expected_path: Path, rtol: float = RTOL, atol: float = ATOL) -> None:
    actual_rows = _read_csv_rows(actual_path)
    expected_rows = _read_csv_rows(expected_path)
    assert len(actual_rows) == len(expected_rows), (
        f"row count mismatch: {len(actual_rows)} actual vs {len(expected_rows)} expected "
        f"({actual_path} vs {expected_path})"
    )
    assert actual_rows[0].keys() == expected_rows[0].keys(), "CSV column mismatch"

    for i, (a_row, e_row) in enumerate(zip(actual_rows, expected_rows)):
        for col in e_row:
            a_val, e_val = a_row[col], e_row[col]
            if col in _EXACT_CSV_COLUMNS:
                assert a_val == e_val, f"row {i} column {col!r}: {a_val!r} != {e_val!r}"
            else:
                a_f, e_f = float(a_val), float(e_val)
                assert np.allclose(a_f, e_f, rtol=rtol, atol=atol), (
                    f"row {i} column {col!r}: {a_f!r} vs {e_f!r} (rtol={rtol}, atol={atol})"
                )


def compare_npz_state(actual_path: Path, expected_path: Path, rtol: float = RTOL, atol: float = ATOL) -> None:
    a = np.load(actual_path, allow_pickle=True)
    e = np.load(expected_path, allow_pickle=True)
    for key in ("xs", "ys", "theta", "omega"):
        assert key in e.files, f"reference NPZ missing key {key!r}"
        assert key in a.files, f"actual NPZ missing key {key!r}"
        assert np.allclose(a[key], e[key], rtol=rtol, atol=atol), f"NPZ array {key!r} differs beyond tolerance"


def compare_metadata_derived(actual_path: Path, expected_path: Path, rtol: float = RTOL, atol: float = ATOL) -> None:
    """Compare the physically meaningful scalar fields of metadata JSON.

    Deliberately excludes fields that are timestamps, wall-clock runtime, or
    absolute file paths (those legitimately differ every run).
    """
    with open(actual_path) as fh:
        a_meta = json.load(fh)
    with open(expected_path) as fh:
        e_meta = json.load(fh)

    a_derived, e_derived = a_meta["derived"], e_meta["derived"]
    assert a_derived.keys() == e_derived.keys(), "metadata 'derived' key set changed"
    for key, e_val in e_derived.items():
        a_val = a_derived[key]
        if isinstance(e_val, (int, float)) and not isinstance(e_val, bool):
            assert np.allclose(a_val, e_val, rtol=rtol, atol=atol), (
                f"derived[{key!r}]: {a_val!r} vs {e_val!r}"
            )
        else:
            assert a_val == e_val, f"derived[{key!r}]: {a_val!r} != {e_val!r}"

    assert a_meta["config"] == e_meta["config"], "metadata 'config' block changed"
