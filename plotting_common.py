"""plotting_common.py — shared needle-lattice drawing code.

Previously ~150 lines (needle_halves, color/style constants, the
polygon+pivot-circle drawing loop) were duplicated verbatim between
compass_generate_images.py and plot_default_geometries_clean.py (see
docs/AUDIT.md P2 item 9). Both now import from here instead.
"""

from __future__ import annotations

import numpy as np
from matplotlib.patches import Circle, Polygon

from sim_config import load_config

_CFG = load_config()
_RENDER_PHYS = _CFG.physics.needle_render
_RENDER_NUM = _CFG.numerics.rendering

BLUE_NORTH = _RENDER_PHYS.colors["blue_north"]
WHITE_SOUTH = _RENDER_PHYS.colors["white_south"]
EDGE = _RENDER_PHYS.colors["edge"]
PIVOT = _RENDER_PHYS.colors["pivot"]
BACKGROUND = _RENDER_PHYS.colors["background"]

PIVOT_RADIUS_FRAC = _RENDER_PHYS.pivot_radius_frac
PIVOT_INNER_RADIUS_FRAC = _RENDER_PHYS.pivot_inner_radius_frac

DEFAULT_DPI = _RENDER_NUM.dpi_default
DEFAULT_FIGSIZE = _RENDER_NUM.figsize_default


def needle_halves(x, y, theta, length, width):
    """Return blue/north and white/south triangular halves of one needle.

    For theta = 0: blue/north side points left, white/south side points right.
    """
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


def draw_lattice(ax, xs, ys, theta, needle_len, needle_width, r_nn, *, clean=True):
    """Draw every needle (two colored triangular halves + a pivot circle)
    into `ax`, then set its limits/aspect.

    `clean=True` hides the axes entirely; `clean=False` shows axis
    labels/grid instead (used by callers that want a labeled figure).
    """
    for x, y, th in zip(xs, ys, theta):
        north_half, south_half = needle_halves(x, y, th, needle_len, needle_width)

        ax.add_patch(
            Polygon(north_half, closed=True, facecolor=BLUE_NORTH, edgecolor=EDGE,
                    linewidth=0.8, joinstyle="miter", zorder=2)
        )
        ax.add_patch(
            Polygon(south_half, closed=True, facecolor=WHITE_SOUTH, edgecolor=EDGE,
                    linewidth=0.8, joinstyle="miter", zorder=2)
        )
        ax.add_patch(
            Circle((x, y), PIVOT_RADIUS_FRAC * r_nn, facecolor=PIVOT,
                   edgecolor=EDGE, linewidth=0.6, zorder=5)
        )
        ax.add_patch(
            Circle((x, y), PIVOT_INNER_RADIUS_FRAC * r_nn, facecolor=PIVOT,
                   edgecolor="white", linewidth=0.35, zorder=6)
        )

    margin = 0.8 * r_nn
    ax.set_xlim(xs.min() - margin, xs.max() + margin)
    ax.set_ylim(ys.min() - margin, ys.max() + margin)
    ax.set_aspect("equal")

    if clean:
        ax.axis("off")
    else:
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.grid(alpha=0.25)
