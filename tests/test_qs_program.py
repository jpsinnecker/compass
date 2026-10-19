"""Quasi-static field program (--qs_program) and return-point-memory metrics (rpm.py)."""
import csv
import json

import numpy as np
import pytest

import compass
import rpm


def _run(tmp_path, tag, *extra):
    argv = ["--out_dir", str(tmp_path), "--tag", tag, "--png_dpi", "20", "--seed", "11",
            "--damping_noise", "0", "--field_mode", "hysteresis_qs", "--qs_dB_factor", "0.05",
            "--qs_tol_torque", "1e-4", "--qs_tol_omega", "1e-4", *extra]
    compass.main(argv)
    with open(tmp_path / "data" / f"{tag}.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    with open(tmp_path / "meta" / f"{tag}.json") as fh:
        meta = json.load(fh)
    return rows, meta


def test_program_path_legs_and_exact_turning_points():
    path = compass.build_qs_program_path([1.0, -0.3, 0.2, -0.3], 0.1)
    labels = [b for _, b in path]
    assert labels[0] == "L01_up" and labels[-1] == "L04_down"
    assert [b for b in dict.fromkeys(labels)] == ["L01_up", "L02_down", "L03_up", "L04_down"]
    ends = [path[i][0] for i in range(len(path)) if i == len(path) - 1 or path[i + 1][1] != path[i][1]]
    assert ends == [1.0, -0.3, 0.2, -0.3]               # turning points reached exactly
    steps = np.abs(np.diff([0.0] + [B for B, _ in path]))
    assert np.all(steps <= 0.1 + 1e-12) and np.all(steps > 0.05)
    with pytest.raises(ValueError):
        compass.build_qs_program_path([1.0, 1.0], 0.1)
    assert compass.parse_qs_program("4, -0.3; 0.2") == [4.0, -0.3, 0.2]


def test_program_reproduces_standard_loop(tmp_path):
    # The five-branch loop written as a program gives the same rest states.
    common = ["--geometry", "square", "--N", "3", "--M", "3", "--B_max_factor", "2.5"]
    rows_std, _ = _run(tmp_path, "std", *common)
    rows_prg, meta = _run(tmp_path, "prg", *common, "--qs_program", "2.5,0,-2.5,0,2.5")
    assert len(rows_std) == len(rows_prg)
    m_std = np.array([float(r["M_proj"]) for r in rows_std])
    m_prg = np.array([float(r["M_proj"]) for r in rows_prg])
    assert np.max(np.abs(m_std - m_prg)) < 1e-6
    assert meta["quasi_static"]["program_Bref"] == [2.5, 0.0, -2.5, 0.0, 2.5]
    assert meta["energy_ledger"]["balance_residual_rel"] < 1e-3


def test_program_needs_quasi_static_mode(tmp_path):
    with pytest.raises(SystemExit):
        compass.main(["--out_dir", str(tmp_path), "--field_mode", "hysteresis", "--qs_program", "1,-1"])


def test_rpm_metrics_on_minor_loop(tmp_path):
    # Minor excursion that reverses nothing: the return state equals the first one.
    _run(tmp_path, "minor", "--geometry", "square", "--N", "3", "--M", "3", "--noise", "0.3",
         "--qs_program", "2.5,-0.4,0.2,-0.4", "--qs_state_log")
    run = rpm.load_run(str(tmp_path), "minor")
    tp = rpm.turning_points(run)
    assert list(tp.leg) == [1, 2, 3, 4]
    assert list(np.round(tp.B_Bref, 9)) == [2.5, -0.4, 0.2, -0.4]
    ret = rpm.rpm_returns(run)
    assert len(ret) == 1 and ret.first_leg[0] == 2 and ret.return_leg[0] == 4
    # Rest tolerance 1e-4 here (1e-5 in production): angles agree to ~1e-3 rad.
    assert ret.flip_mismatch[0] == 0.0 and ret.rms_rad[0] < 1e-2


def test_state_distance():
    th = np.array([0.1, 3.0, -3.1])
    d = rpm.state_distance(th, th + 2 * np.pi, 0.0)
    assert d["rms_rad"] < 1e-12 and d["flip_mismatch"] == 0.0
    d = rpm.state_distance(np.array([0.0, 0.0]), np.array([0.0, np.pi]), 0.0)
    assert d["flip_mismatch"] == 0.5 and d["dM"] == pytest.approx(1.0)
