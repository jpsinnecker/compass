#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from sim_config import load_config
from plotting_common import BACKGROUND, DEFAULT_DPI, DEFAULT_FIGSIZE, draw_lattice

_CFG = load_config()


def load_state(npz_path):
    """Load (xs, ys, theta, needle_len, needle_width, r_nn) from a run's
    final-state NPZ. Tolerant of older/incomplete metadata (see
    docs/AUDIT.md bug B5): falls back to the same config-driven geometry
    guesses compass_generate_images.py's load_npz_state() uses, instead of
    raising KeyError when metadata_json or a 'derived' field is missing.
    """
    npz_path = Path(npz_path)

    if not npz_path.exists():
        raise FileNotFoundError(f"Could not find: {npz_path}")

    data = np.load(npz_path, allow_pickle=True)

    xs = data["xs"]
    ys = data["ys"]
    theta = data["theta"]

    metadata = {}
    if "metadata_json" in data.files:
        metadata = json.loads(str(data["metadata_json"]))
    derived = metadata.get("derived", {})

    # Width/length ratio guess uses compass_engine's legacy_width_to_length_ratio
    # (0.22), matching --use_legacy_size_from_R -- not the unrelated 0.30
    # research default (see docs/AUDIT.md P2 item 10).
    _fallback_cfg = _CFG.physics.compass_generate_images
    r_nn = float(data["r_nn"]) if "r_nn" in data.files else float(derived.get("r_nn_m", _fallback_cfg.r_nn_fallback))
    needle_len = float(derived.get("needle_len_m", _fallback_cfg.needle_len_to_r_nn_fallback_ratio * r_nn))
    needle_width = float(derived.get("needle_width_m", _CFG.physics.compass_engine.legacy_width_to_length_ratio * needle_len))

    return xs, ys, theta, needle_len, needle_width, r_nn


def plot_geometry(npz_path, out_png, dpi=DEFAULT_DPI,
                  figsize=DEFAULT_FIGSIZE, transparent=False):
    xs, ys, theta, needle_len, needle_width, r_nn = load_state(npz_path)

    fig, ax = plt.subplots(
        figsize=(figsize, figsize),
        facecolor="none" if transparent else BACKGROUND,
    )

    draw_lattice(ax, xs, ys, theta, needle_len, needle_width, r_nn, clean=True)

    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        out_png,
        dpi=dpi,
        bbox_inches="tight",
        pad_inches=0,
        facecolor="none" if transparent else BACKGROUND,
        transparent=transparent,
    )
    plt.close(fig)

    print(f"saved {out_png}")


def find_single_final_npz(input_dir):
    states_dir = Path(input_dir) / "states"

    if not states_dir.exists():
        raise FileNotFoundError(f"Could not find states directory: {states_dir}")

    candidates = sorted(states_dir.glob("*_final.npz"))

    if not candidates:
        raise FileNotFoundError(f"No *_final.npz files found in: {states_dir}")

    if len(candidates) > 1:
        print("Multiple *_final.npz files found. Using the first one:")
        for c in candidates:
            print(f"  {c}")

    return candidates[0]


def find_default_geometry_files(input_dir, tag_suffix="_default"):
    input_dir = Path(input_dir)
    mapping = {}

    for geometry in ["square", "triangular", "honeycomb"]:
        expected = (
            input_dir
            / geometry
            / "states"
            / f"{geometry}{tag_suffix}_final.npz"
        )

        if expected.exists():
            mapping[geometry] = expected
            continue

        candidates = sorted((input_dir / geometry / "states").glob("*_final.npz"))

        if candidates:
            mapping[geometry] = candidates[0]
            continue

        raise FileNotFoundError(
            f"Could not find final state for {geometry}.\n"
            f"Expected: {expected}"
        )

    return mapping


def main():
    parser = argparse.ArgumentParser(
        description="Generate clean PNG geometry figures from compass.py .npz output."
    )
    run_cfg = _CFG.run.plot_default_geometries_clean

    parser.add_argument(
        "--input_dir",
        type=str,
        default=run_cfg.input_dir,
        help="Input directory. Default: geometry_figures",
    )

    parser.add_argument(
        "--input_file",
        type=str,
        default=run_cfg.input_file,
        help="Specific .npz final-state file to plot.",
    )

    parser.add_argument(
        "--output_file",
        type=str,
        default=run_cfg.output_file,
        help="Output PNG file for --input_file or --single.",
    )

    parser.add_argument(
        "--output_dir",
        type=str,
        default=run_cfg.output_dir,
        help="Output directory for batch mode.",
    )

    parser.add_argument(
        "--single",
        action="store_true",
        default=run_cfg.single,
        help="Treat --input_dir as one compass.py output directory.",
    )

    parser.add_argument(
        "--tag_suffix",
        type=str,
        default=run_cfg.tag_suffix,
        help="Suffix used for default batch files. Default: _default",
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=DEFAULT_DPI,
        help=f"PNG resolution. Default: {DEFAULT_DPI}",
    )

    parser.add_argument(
        "--figsize",
        type=float,
        default=DEFAULT_FIGSIZE,
        help=f"Square figure size in inches. Default: {DEFAULT_FIGSIZE}",
    )

    parser.add_argument(
        "--transparent",
        action="store_true",
        default=run_cfg.transparent,
        help="Save PNG with transparent background.",
    )

    args = parser.parse_args()

    if args.input_file is not None:
        input_file = Path(args.input_file)
        output_file = (
            Path(args.output_file)
            if args.output_file is not None
            else input_file.with_name(input_file.stem + "_clean.png")
        )

        plot_geometry(
            input_file,
            output_file,
            dpi=args.dpi,
            figsize=args.figsize,
            transparent=args.transparent,
        )
        return

    if args.single:
        input_dir = Path(args.input_dir)
        input_file = find_single_final_npz(input_dir)
        output_file = (
            Path(args.output_file)
            if args.output_file is not None
            else input_dir / (input_file.stem.replace("_final", "") + "_clean.png")
        )

        plot_geometry(
            input_file,
            output_file,
            dpi=args.dpi,
            figsize=args.figsize,
            transparent=args.transparent,
        )
        return

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir) if args.output_dir else input_dir

    mapping = find_default_geometry_files(input_dir, tag_suffix=args.tag_suffix)

    for geometry, input_file in mapping.items():
        output_file = output_dir / f"{geometry}_default_geometry_clean_colored.png"
        plot_geometry(
            input_file,
            output_file,
            dpi=args.dpi,
            figsize=args.figsize,
            transparent=args.transparent,
        )


if __name__ == "__main__":
    main()
