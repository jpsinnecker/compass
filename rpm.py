#!/usr/bin/env python3
"""
rpm.py -- return-point-memory (RPM) metrics for quasi-static field programs
(compass.py --field_mode hysteresis_qs --qs_program ... --qs_state_log).

A run visits turning points in order (legs Lkk_up / Lkk_down). Return-point
memory means that when the field comes back to an earlier turning point after
a minor excursion, the array is in the same state as at the first passage, and
the outer branch then continues as if the excursion had not happened.

Metrics between two states theta1, theta2 (radians, one value per needle):
  rms_rad        rms of the wrapped angle difference
  flip_mismatch  fraction of needles whose projection on the drive axis has a
                 different sign (the needle is reversed in one state only)
  dM             difference of the projected magnetization per needle

Usage as a script:
  python3 rpm.py <run_dir> [<run_dir> ...]       turning-point table per run
  python3 rpm.py <run_dir> --ref <ref_run_dir> --legs 5:2
                 compare leg 5 of the run with leg 2 of a reference run at the
                 same field values (outer branch after the excursion)
"""
from __future__ import annotations

import argparse
import glob
import json
import os
from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class QSRun:
    tag: str
    df: pd.DataFrame          # one row per field step (row 0 = initial state)
    theta: np.ndarray         # rest-state angles [row, needle]
    phi_ext: float            # drive-axis angle [rad]
    meta: dict


def load_run(run_dir: str, tag: str | None = None) -> QSRun:
    """Load CSV, metadata and the per-step state log of one quasi-static run."""
    metas = sorted(glob.glob(os.path.join(run_dir, "meta", "*.json")))
    if tag is None:
        if len(metas) != 1:
            raise ValueError(f"{run_dir}: {len(metas)} runs, pass tag")
        tag = os.path.splitext(os.path.basename(metas[0]))[0]
    with open(os.path.join(run_dir, "meta", f"{tag}.json")) as fh:
        meta = json.load(fh)
    df = pd.read_csv(os.path.join(run_dir, "data", f"{tag}.csv"))
    z = np.load(os.path.join(run_dir, "states", f"{tag}_qs_states.npz"))
    theta = np.asarray(z["theta"])
    if theta.shape[0] != len(df):
        raise ValueError(f"{tag}: state log has {theta.shape[0]} rows, CSV {len(df)}")
    return QSRun(tag, df, theta, float(np.asarray(z["phi_ext"]).ravel()[0]), meta)


def state_distance(th1: np.ndarray, th2: np.ndarray, phi_ext: float) -> dict:
    d = np.angle(np.exp(1j * (np.asarray(th1) - np.asarray(th2))))
    p1 = np.cos(np.asarray(th1) - phi_ext)
    p2 = np.cos(np.asarray(th2) - phi_ext)
    return dict(rms_rad=float(np.sqrt(np.mean(d**2))),
                flip_mismatch=float(np.mean(np.sign(p1) != np.sign(p2))),
                dM=float(np.mean(p1) - np.mean(p2)))


def turning_points(run: QSRun) -> pd.DataFrame:
    """Row index, leg and field (units of B_ref) of the last row of every leg."""
    br = run.df.branch.astype(str).values
    B_ref = run.meta["derived"]["B_ref_T"]
    rows = []
    for i in range(1, len(br)):
        if br[i].startswith("L") and (i == len(br) - 1 or br[i + 1] != br[i]):
            rows.append(dict(row=i, leg=int(br[i][1:3]), B_Bref=float(run.df.B_scalar_T.iloc[i]) / B_ref))
    return pd.DataFrame(rows)


def rpm_returns(run: QSRun, tol: float = 1e-9) -> pd.DataFrame:
    """Compare every turning point with the first earlier turning point at the same field."""
    tp = turning_points(run)
    out = []
    for j in range(len(tp)):
        for i in range(j):
            if abs(tp.B_Bref[i] - tp.B_Bref[j]) <= tol * max(1.0, abs(tp.B_Bref[i])):
                d = state_distance(run.theta[tp.row[j]], run.theta[tp.row[i]], run.phi_ext)
                out.append(dict(B_Bref=tp.B_Bref[j], first_leg=int(tp.leg[i]), return_leg=int(tp.leg[j]), **d))
                break
    return pd.DataFrame(out)


def compare_legs(run: QSRun, leg: int, ref: QSRun, ref_leg: int, tol_frac: float = 0.1) -> pd.DataFrame:
    """State distance along a leg of `run` versus a leg of `ref`, at equal field values.

    Fields are matched within tol_frac of the step; legs on a common grid
    (turning points multiples of dB) match point by point.
    """
    B_ref = run.meta["derived"]["B_ref_T"]
    dB = run.meta["quasi_static"]["dB_T"]
    a = run.df.index[run.df.branch == f"L{leg:02d}_{'up' if _leg_dir(run, leg) > 0 else 'down'}"]
    b = ref.df.index[ref.df.branch.astype(str).str.startswith(f"L{ref_leg:02d}_")]
    Bb = ref.df.B_scalar_T.values[b]
    out = []
    for i in a:
        Bi = run.df.B_scalar_T.values[i]
        k = int(np.argmin(np.abs(Bb - Bi)))
        if abs(Bb[k] - Bi) <= tol_frac * dB:
            d = state_distance(run.theta[i], ref.theta[b[k]], run.phi_ext)
            out.append(dict(B_Bref=Bi / B_ref, **d))
    return pd.DataFrame(out)


def _leg_dir(run: QSRun, leg: int) -> int:
    lab = run.df.branch.astype(str)
    hit = lab[lab.str.startswith(f"L{leg:02d}_")]
    if hit.empty:
        raise ValueError(f"{run.tag}: no leg {leg}")
    return 1 if hit.iloc[0].endswith("_up") else -1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dirs", nargs="+")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--ref", default=None, help="reference run directory for --legs")
    ap.add_argument("--legs", default=None, help="LEG:REF_LEG, e.g. 5:2")
    args = ap.parse_args(argv)
    pd.set_option("display.width", 160)
    for rd in args.run_dirs:
        run = load_run(rd, args.tag)
        print(f"== {run.tag}  program {run.meta['quasi_static'].get('program_Bref')}")
        print(turning_points(run).to_string(index=False))
        ret = rpm_returns(run)
        print(ret.to_string(index=False) if not ret.empty else "(no returns to an earlier turning point)")
        if args.ref and args.legs:
            leg, ref_leg = (int(v) for v in args.legs.split(":"))
            cmp = compare_legs(run, leg, load_run(args.ref), ref_leg)
            print(f"-- leg {leg} vs reference leg {ref_leg}: max flip_mismatch {cmp.flip_mismatch.max():.4f}, "
                  f"max |dM| {cmp.dM.abs().max():.4f} over {len(cmp)} fields")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
