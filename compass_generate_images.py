#!/usr/bin/env python3
"""
compass_generate_images.py

Post-process one or more compass.py run directories and regenerate all standard
PNG images from the saved CSV and NPZ files.

The script does not import compass.py. It only reads:

    <run_dir>/states/<tag>_initial.npz
    <run_dir>/states/<tag>_final.npz
    <run_dir>/data/<tag>.csv
    <run_dir>/meta/<tag>.json        optional

Images written:

    <run_dir>/images/<tag>_final_lattice.png
    <run_dir>/images/<tag>_initial_final.png
    <run_dir>/images/<tag>_quicklook.png
    <run_dir>/images/<tag>_field_protocol.png
    <run_dir>/images/<tag>_energies.png
    <run_dir>/images/<tag>_order_parameters.png

Examples
--------

Process one run directory:

    python3 compass_generate_images.py --run_dir test_honeycomb

Process all run directories below field_mode_tests:

    python3 compass_generate_images.py --run_dir field_mode_tests --recursive

Use transparent background for lattice images:

    python3 compass_generate_images.py --run_dir test_honeycomb --transparent

Include axes in lattice images:

    python3 compass_generate_images.py --run_dir test_honeycomb --with_axes
"""

import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from sim_config import load_config
from plotting_common import BACKGROUND, DEFAULT_DPI, draw_lattice

_CFG = load_config()


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

def ensure_dir(path):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_csv_structured(csv_path):
    """
    Read compass.py CSV output as a structured NumPy array.
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(csv_path)

    data = np.genfromtxt(
        csv_path,
        delimiter=",",
        names=True,
        dtype=None,
        encoding=None,
        autostrip=True,
    )

    # If the CSV contains one data row, genfromtxt returns a scalar structured
    # array. Convert it to a one-element array for uniform indexing.
    if data.shape == ():
        data = np.array([data], dtype=data.dtype)

    return data


def col(data, name, default=None):
    """
    Return one named CSV column as a NumPy array, or a default array.
    """
    if name in data.dtype.names:
        return np.asarray(data[name])

    if default is None:
        raise KeyError(f"CSV column not found: {name}")

    n = len(data)
    return np.full(n, default)


def numeric_col(data, name, default=np.nan):
    """
    Return one named CSV column converted to float.
    """
    return np.asarray(col(data, name, default), dtype=float)


def load_npz_state(npz_path):
    """
    Load one compass.py state NPZ file.
    """
    npz_path = Path(npz_path)
    if not npz_path.exists():
        raise FileNotFoundError(npz_path)

    d = np.load(npz_path, allow_pickle=True)

    xs = np.asarray(d["xs"], dtype=float)
    ys = np.asarray(d["ys"], dtype=float)
    theta = np.asarray(d["theta"], dtype=float)
    omega = np.asarray(d["omega"], dtype=float) if "omega" in d.files else np.zeros_like(theta)

    metadata = {}
    if "metadata_json" in d.files:
        metadata = json.loads(str(d["metadata_json"]))

    derived = metadata.get("derived", {})
    config = metadata.get("config", {})

    # Fallbacks keep the script usable with older NPZ files. The width/length
    # ratio guess uses compass_engine's legacy_width_to_length_ratio (0.22),
    # matching compass.py's --use_legacy_size_from_R convention that these
    # older files most likely used -- not the unrelated 0.30 research default
    # (see docs/AUDIT.md P2 item 10).
    _fallback_cfg = _CFG.physics.compass_generate_images
    r_nn = float(d["r_nn"]) if "r_nn" in d.files else float(derived.get("r_nn_m", _fallback_cfg.r_nn_fallback))
    needle_len = float(derived.get("needle_len_m", _fallback_cfg.needle_len_to_r_nn_fallback_ratio * r_nn))
    needle_width = float(derived.get("needle_width_m", _CFG.physics.compass_engine.legacy_width_to_length_ratio * needle_len))

    return {
        "xs": xs,
        "ys": ys,
        "theta": theta,
        "omega": omega,
        "r_nn": r_nn,
        "needle_len": needle_len,
        "needle_width": needle_width,
        "metadata": metadata,
        "geometry": config.get("geometry", npz_path.stem),
        "field_mode": config.get("field_mode", ""),
    }


def find_run_dirs(root, recursive=False):
    """
    Find directories that look like compass.py output directories.
    A run directory has either states/ or data/.
    """
    root = Path(root)

    if not recursive:
        return [root]

    dirs = []
    for p in [root] + list(root.rglob("*")):
        if not p.is_dir():
            continue
        if (p / "states").exists() or (p / "data").exists():
            dirs.append(p)

    # Remove nested states/data/images/meta directories themselves.
    skip_names = {"states", "data", "images", "meta"}
    dirs = [d for d in dirs if d.name not in skip_names]

    # Deduplicate while preserving order.
    out = []
    seen = set()
    for d in dirs:
        rp = d.resolve()
        if rp not in seen:
            seen.add(rp)
            out.append(d)

    return out


def discover_tags(run_dir):
    """
    Discover tags from states/*_final.npz and data/*.csv.
    """
    run_dir = Path(run_dir)
    tags = set()

    states_dir = run_dir / "states"
    data_dir = run_dir / "data"

    if states_dir.exists():
        for p in states_dir.glob("*_final.npz"):
            tags.add(p.name.replace("_final.npz", ""))

    if data_dir.exists():
        for p in data_dir.glob("*.csv"):
            tags.add(p.stem)

    return sorted(tags)


# ---------------------------------------------------------------------------
# Lattice drawing
# ---------------------------------------------------------------------------

def draw_lattice_state(ax, state, theta=None, clean=True):
    """Draw a run's lattice `state` dict (see load_npz_state) into `ax`."""
    theta_values = state["theta"] if theta is None else theta
    draw_lattice(
        ax, state["xs"], state["ys"], theta_values,
        state["needle_len"], state["needle_width"], state["r_nn"],
        clean=clean,
    )


def plot_final_lattice(final_npz, out_png, *, dpi=DEFAULT_DPI,
                       transparent=False, clean=True):
    state = load_npz_state(final_npz)

    fig, ax = plt.subplots(
        figsize=(6, 6),
        facecolor="none" if transparent else BACKGROUND,
    )

    draw_lattice_state(ax, state, clean=clean)

    if not clean:
        title = f"{state['geometry']} final"
        if state["field_mode"]:
            title += f" — {state['field_mode']}"
        ax.set_title(title)

    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    out_png = Path(out_png)
    ensure_dir(out_png.parent)
    fig.savefig(
        out_png,
        dpi=dpi,
        bbox_inches="tight",
        pad_inches=0 if clean else 0.1,
        facecolor="none" if transparent else BACKGROUND,
        transparent=transparent,
    )
    plt.close(fig)
    print(f"saved {out_png}")


def plot_initial_final(initial_npz, final_npz, out_png, *, dpi=DEFAULT_DPI,
                       transparent=False, clean=True, panel_titles=True):
    s0 = load_npz_state(initial_npz)
    s1 = load_npz_state(final_npz)

    # Use final metadata/geometry, but initial theta.
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(10, 5),
        facecolor="none" if transparent else BACKGROUND,
    )

    draw_lattice_state(axes[0], s1, theta=s0["theta"], clean=clean)
    draw_lattice_state(axes[1], s1, theta=s1["theta"], clean=clean)

    if panel_titles:
        axes[0].set_title("Initial state")
        axes[1].set_title("Final state")

    if clean:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0.02)
    else:
        fig.tight_layout()

    out_png = Path(out_png)
    ensure_dir(out_png.parent)
    fig.savefig(
        out_png,
        dpi=dpi,
        bbox_inches="tight",
        pad_inches=0 if clean else 0.1,
        facecolor="none" if transparent else BACKGROUND,
        transparent=transparent,
    )
    plt.close(fig)
    print(f"saved {out_png}")


# ---------------------------------------------------------------------------
# CSV plots
# ---------------------------------------------------------------------------

def plot_quicklook(csv_path, out_png, *, dpi=DEFAULT_DPI):
    data = read_csv_structured(csv_path)

    t = numeric_col(data, "t_s")
    B = numeric_col(data, "B_scalar_T")
    M = numeric_col(data, "M_proj")
    S1 = numeric_col(data, "S1")
    S2 = numeric_col(data, "S2")
    Etot = numeric_col(data, "E_total_J")
    omega = numeric_col(data, "omega_rms_rad_s")

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))

    ax = axes[0, 0]
    if np.nanmax(B) > np.nanmin(B):
        ax.plot(B * 1e3, M, marker=".", linewidth=1)
        ax.set_xlabel("B scalar (mT)")
    else:
        ax.plot(t, M, marker=".", linewidth=1)
        ax.set_xlabel("t (s)")
    ax.set_ylabel("M_proj")
    ax.set_title("Magnetization")

    ax = axes[0, 1]
    ax.plot(t, S1, label="S1", linewidth=1)
    ax.plot(t, S2, label="S2", linewidth=1)
    ax.set_xlabel("t (s)")
    ax.set_ylabel("order")
    ax.set_title("Order parameters")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.plot(t, Etot, linewidth=1)
    ax.set_xlabel("t (s)")
    ax.set_ylabel("E_total (J)")
    ax.set_title("Total energy")

    ax = axes[1, 1]
    ax.plot(t, omega, linewidth=1)
    ax.set_xlabel("t (s)")
    ax.set_ylabel("omega RMS (rad/s)")
    ax.set_title("Angular velocity")

    for ax in axes.flat:
        ax.grid(alpha=0.25)

    fig.tight_layout()
    out_png = Path(out_png)
    ensure_dir(out_png.parent)
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)
    print(f"saved {out_png}")


def plot_field_protocol(csv_path, out_png, *, dpi=DEFAULT_DPI):
    data = read_csv_structured(csv_path)

    t = numeric_col(data, "t_s")
    Bx = numeric_col(data, "Bx_T")
    By = numeric_col(data, "By_T")
    Bs = numeric_col(data, "B_scalar_T")
    M = numeric_col(data, "M_proj")

    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    axes[0].plot(t, Bx * 1e3, label="Bx", linewidth=1)
    axes[0].plot(t, By * 1e3, label="By", linewidth=1)
    axes[0].plot(t, Bs * 1e3, label="B scalar", linewidth=1, linestyle="--")
    axes[0].set_ylabel("Field (mT)")
    axes[0].set_title("Field protocol")
    axes[0].legend(frameon=False)

    axes[1].plot(t, M, linewidth=1)
    axes[1].set_xlabel("t (s)")
    axes[1].set_ylabel("M_proj")
    axes[1].set_title("Projected magnetization")

    for ax in axes:
        ax.grid(alpha=0.25)

    fig.tight_layout()
    out_png = Path(out_png)
    ensure_dir(out_png.parent)
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)
    print(f"saved {out_png}")


def plot_energies(csv_path, out_png, *, dpi=DEFAULT_DPI):
    data = read_csv_structured(csv_path)

    t = numeric_col(data, "t_s")
    Edip = numeric_col(data, "E_dip_J")
    Eext = numeric_col(data, "E_ext_J")
    Ekin = numeric_col(data, "E_kin_J")
    Etot = numeric_col(data, "E_total_J")

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(t, Edip, label="E_dip", linewidth=1)
    ax.plot(t, Eext, label="E_ext", linewidth=1)
    ax.plot(t, Ekin, label="E_kin", linewidth=1)
    ax.plot(t, Etot, label="E_total", linewidth=1.3)

    ax.set_xlabel("t (s)")
    ax.set_ylabel("Energy (J)")
    ax.set_title("Energy terms")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)

    fig.tight_layout()
    out_png = Path(out_png)
    ensure_dir(out_png.parent)
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)
    print(f"saved {out_png}")


def plot_order_parameters(csv_path, out_png, *, dpi=DEFAULT_DPI):
    data = read_csv_structured(csv_path)

    t = numeric_col(data, "t_s")
    S1 = numeric_col(data, "S1")
    S2 = numeric_col(data, "S2")
    q = numeric_col(data, "q_axis")
    ff = numeric_col(data, "flip_field", default=0.0)
    fa = numeric_col(data, "flip_angle", default=0.0)

    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    axes[0].plot(t, S1, label="S1", linewidth=1)
    axes[0].plot(t, S2, label="S2", linewidth=1)
    axes[0].plot(t, q, label="q_axis", linewidth=1)
    axes[0].set_ylabel("Order")
    axes[0].set_title("Order parameters")
    axes[0].legend(frameon=False)

    axes[1].plot(t, ff, label="flip_field", linewidth=1)
    axes[1].plot(t, fa, label="flip_angle", linewidth=1)
    axes[1].set_xlabel("t (s)")
    axes[1].set_ylabel("Flip count per log interval")
    axes[1].set_title("Avalanche proxies")
    axes[1].legend(frameon=False)

    for ax in axes:
        ax.grid(alpha=0.25)

    fig.tight_layout()
    out_png = Path(out_png)
    ensure_dir(out_png.parent)
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)
    print(f"saved {out_png}")


# ---------------------------------------------------------------------------
# Run processing
# ---------------------------------------------------------------------------

def process_tag(run_dir, tag, args):
    run_dir = Path(run_dir)
    image_dir = ensure_dir(run_dir / "images")

    final_npz = run_dir / "states" / f"{tag}_final.npz"
    initial_npz = run_dir / "states" / f"{tag}_initial.npz"
    csv_path = run_dir / "data" / f"{tag}.csv"

    if final_npz.exists():
        plot_final_lattice(
            final_npz,
            image_dir / f"{tag}_final_lattice.png",
            dpi=args.dpi,
            transparent=args.transparent,
            clean=not args.with_axes,
        )

    if initial_npz.exists() and final_npz.exists():
        plot_initial_final(
            initial_npz,
            final_npz,
            image_dir / f"{tag}_initial_final.png",
            dpi=args.dpi,
            transparent=args.transparent,
            clean=not args.with_axes,
            panel_titles=not args.no_panel_titles,
        )

    if csv_path.exists():
        plot_quicklook(
            csv_path,
            image_dir / f"{tag}_quicklook.png",
            dpi=args.dpi,
        )
        plot_field_protocol(
            csv_path,
            image_dir / f"{tag}_field_protocol.png",
            dpi=args.dpi,
        )
        plot_energies(
            csv_path,
            image_dir / f"{tag}_energies.png",
            dpi=args.dpi,
        )
        plot_order_parameters(
            csv_path,
            image_dir / f"{tag}_order_parameters.png",
            dpi=args.dpi,
        )


def process_run_dir(run_dir, args):
    run_dir = Path(run_dir)
    tags = discover_tags(run_dir)

    if not tags:
        if args.verbose:
            print(f"[skip] no compass.py outputs found in {run_dir}")
        return

    if args.verbose:
        print(f"[run] {run_dir}")
        print(f"      tags: {', '.join(tags)}")

    for tag in tags:
        process_tag(run_dir, tag, args)


def build_parser():
    p = argparse.ArgumentParser(
        description="Generate PNG images from compass.py CSV and NPZ outputs."
    )
    run_cfg = _CFG.run.compass_generate_images

    p.add_argument(
        "--run_dir",
        type=str,
        default=run_cfg.run_dir,
        help="Compass.py output directory. Default: current directory.",
    )

    p.add_argument(
        "--recursive",
        action="store_true",
        default=run_cfg.recursive,
        help="Process all compass.py run directories below --run_dir.",
    )

    p.add_argument(
        "--dpi",
        type=int,
        default=DEFAULT_DPI,
        help=f"PNG resolution. Default: {DEFAULT_DPI}.",
    )

    p.add_argument(
        "--transparent",
        action="store_true",
        default=run_cfg.transparent,
        help="Use transparent background for lattice PNGs.",
    )

    p.add_argument(
        "--with_axes",
        action="store_true",
        default=run_cfg.with_axes,
        help="Show axes/grid on lattice images.",
    )

    p.add_argument(
        "--no_panel_titles",
        action="store_true",
        default=run_cfg.no_panel_titles,
        help="Remove 'Initial state' and 'Final state' titles.",
    )

    p.add_argument(
        "--verbose",
        action="store_true",
        default=run_cfg.verbose,
        help="Print discovered run directories and tags.",
    )

    return p


def main():
    args = build_parser().parse_args()

    for run_dir in find_run_dirs(args.run_dir, recursive=args.recursive):
        process_run_dir(run_dir, args)


if __name__ == "__main__":
    main()
