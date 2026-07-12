#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle

BLUE_NORTH = "#0017B8"
WHITE_SOUTH = "#F2F2F2"
EDGE = "#4D4D4D"
PIVOT = "#777777"
BACKGROUND = "white"

PIVOT_RADIUS_FRAC = 0.085
PIVOT_INNER_RADIUS_FRAC = 0.025

DEFAULT_DPI = 300
DEFAULT_FIGSIZE = 6.0


def needle_halves(x, y, theta, length, width):
    u = np.array([np.cos(theta), np.sin(theta)])
    v = np.array([-np.sin(theta), np.cos(theta)])

    center = np.array([x, y])
    left_tip = center - 0.5 * length * u
    right_tip = center + 0.5 * length * u
    top_mid = center + 0.5 * width * v
    bottom_mid = center - 0.5 * width * v

    north_half = np.array([left_tip, top_mid, bottom_mid])
    south_half = np.array([right_tip, top_mid, bottom_mid])

    return north_half, south_half


def load_state(npz_path):
    npz_path = Path(npz_path)

    if not npz_path.exists():
        raise FileNotFoundError(f"Could not find: {npz_path}")

    data = np.load(npz_path, allow_pickle=True)

    xs = data["xs"]
    ys = data["ys"]
    theta = data["theta"]

    meta = json.loads(str(data["metadata_json"]))
    derived = meta["derived"]

    needle_len = float(derived["needle_len_m"])
    needle_width = float(derived["needle_width_m"])
    r_nn = float(derived["r_nn_m"])

    return xs, ys, theta, needle_len, needle_width, r_nn


def plot_geometry(npz_path, out_png, dpi=DEFAULT_DPI,
                  figsize=DEFAULT_FIGSIZE, transparent=False):
    xs, ys, theta, needle_len, needle_width, r_nn = load_state(npz_path)

    fig, ax = plt.subplots(
        figsize=(figsize, figsize),
        facecolor="none" if transparent else BACKGROUND,
    )

    for x, y, th in zip(xs, ys, theta):
        north_half, south_half = needle_halves(
            x, y, th, needle_len, needle_width
        )

        ax.add_patch(
            Polygon(
                north_half,
                closed=True,
                facecolor=BLUE_NORTH,
                edgecolor=EDGE,
                linewidth=0.8,
                joinstyle="miter",
                zorder=2,
            )
        )

        ax.add_patch(
            Polygon(
                south_half,
                closed=True,
                facecolor=WHITE_SOUTH,
                edgecolor=EDGE,
                linewidth=0.8,
                joinstyle="miter",
                zorder=2,
            )
        )

        ax.add_patch(
            Circle(
                (x, y),
                PIVOT_RADIUS_FRAC * r_nn,
                facecolor=PIVOT,
                edgecolor=EDGE,
                linewidth=0.6,
                zorder=5,
            )
        )

        ax.add_patch(
            Circle(
                (x, y),
                PIVOT_INNER_RADIUS_FRAC * r_nn,
                facecolor=PIVOT,
                edgecolor="white",
                linewidth=0.35,
                zorder=6,
            )
        )

    ax.set_aspect("equal")

    margin = 0.8 * r_nn
    ax.set_xlim(xs.min() - margin, xs.max() + margin)
    ax.set_ylim(ys.min() - margin, ys.max() + margin)

    ax.axis("off")
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
        description="Generate clean PNG geometry figures from compassV02.py .npz output."
    )

    parser.add_argument(
        "--input_dir",
        type=str,
        default="geometry_figures",
        help="Input directory. Default: geometry_figures",
    )

    parser.add_argument(
        "--input_file",
        type=str,
        default=None,
        help="Specific .npz final-state file to plot.",
    )

    parser.add_argument(
        "--output_file",
        type=str,
        default=None,
        help="Output PNG file for --input_file or --single.",
    )

    parser.add_argument(
        "--output_dir",
        type=str,
        default=None,
        help="Output directory for batch mode.",
    )

    parser.add_argument(
        "--single",
        action="store_true",
        help="Treat --input_dir as one compassV02.py output directory.",
    )

    parser.add_argument(
        "--tag_suffix",
        type=str,
        default="_default",
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
