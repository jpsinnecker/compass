#!/usr/bin/env python3
"""
compassV02.py — research-grade simulation of dipolar compass-needle arrays
Version: V2.0.2, July 2026
File timestamp: 2026-07-10 11:50:24 -03

Purpose
-------
This program simulates two-dimensional arrays of classical magnetic compass
needles. Each needle is represented as a point magnetic dipole fixed at a
lattice site and free to rotate in the sample plane. The equation of motion is
Newton's second law for rotation,

    I theta_ddot = tau_dip + tau_ext - b theta_dot,

where the dipolar field is computed from the full tensor of pair interactions
within a controlled cutoff. The code is intended for research studies of
hysteresis, correlated reversal, damping-controlled avalanches, FORC curves,
step/pulse relaxation, and director memory in classical dipolar arrays.

Major improvements relative to compass.py V79
---------------------------------------------
1. Explicit cutoff convention:
   --cutoff_shells specifies the cutoff as a multiple of r_nn.
   --cutoff_m specifies an absolute SI cutoff and overrides cutoff_shells.
   Internally, the cutoff is always stored in metres and both values are
   written to metadata.

2. Research observables are native outputs:
   Mx, My, M_proj, S1, S2, director angle, square-axis population q_axis,
   flip counters, dipolar/external/kinetic energy, omega_rms and omega_max.
   S1 is polar order. S2 is nematic/director order and is essential for
   compass-lattice studies where opposite directions may cancel magnetically.

3. Avalanche activity is exported directly (V2.1 hardened counters):
   flip_field counts COMMITTED sign reversals relative to the drive axis,
   using a hysteretic (Schmitt-trigger) angular dead band of half-width
   --flip_band_deg around the perpendicular, plus a dwell requirement of
   --flip_dwell_T0 natural periods in the new state. A single underdamped
   needle ringing across the axis therefore contributes exactly one event
   per genuine reversal, not one event per zero crossing.
   flip_angle counts COMMITTED rest-angle displacement events: a needle
   registers an event only when its angle has moved by more than
   --flip_angle_deg away from its last settled reference angle AND it has
   settled there (|omega| below --flip_settle_frac * omega0 for the dwell
   window). The reference angle is then updated. This channel is
   geometry-agnostic and detects partial reversals into frustrated local
   minima that never leave the drive-axis dead band.
   With --event_log, every committed event is written to a per-run CSV
   (step, t, needle_id, channel, theta) for offline spatio-temporal
   avalanche clustering.

3b. Numerical stability monitor:
   The local instantaneous stiffness frequency omega_loc = sqrt(m|B_i|/I)
   is monitored every step from the already-computed total field. Steps
   where max_i(omega_loc) * dt exceeds --dt_guard_alpha are counted and
   reported in metadata and at exit. With --dt_guard_substep, flagged
   steps are re-integrated globally with 4 sub-steps (all needles; a
   partial/local substep would desynchronize the coupled field). Note
   that state-dependent stepping formally breaks symplecticity, so the
   guard is off by default; the monitor itself is always on and free.

4. Energy accounting:
   E_dip = -1/2 sum_i m_i · B_i^dip,
   E_ext = -sum_i m_i · B_ext,
   E_kin = 1/2 I sum_i omega_i^2.
   This enables checks of relaxation pathways, dissipation, and numerical
   stability.

5. Efficient torque calculation:
   Axx, Axy and Ayy pair-interaction tensors are precomputed when memory
   allows. Each time step then requires matrix-vector products instead of
   rebuilding pair distances. CPU/NumPy and GPU/CuPy backends are supported.

6. Field protocols:
   static, hysteresis, forc, sine, step_up, step_down, pulse, demag_rot,
   and demag_linear. The code writes consistent CSV and JSON metadata for all
   protocols.

Scientific caveat
-----------------
This is a continuous-angle inertial dipolar model. It is not a literal model of
nanolithographic artificial spin ice, whose islands have shape anisotropy,
micromagnetic reversal modes, switching-field distributions, and material
specific damping. Its research value is to isolate what follows from dipolar
coupling, geometry, inertia, damping, and field history alone.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import time
from datetime import datetime, timezone
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

try:
    import cupy as cp  # type: ignore
    _CUPY_IMPORT_OK = True
    _CUPY_ERROR = ""
except Exception as exc:  # pragma: no cover - environment dependent
    cp = None  # type: ignore
    _CUPY_IMPORT_OK = False
    _CUPY_ERROR = str(exc).splitlines()[-1] if str(exc) else type(exc).__name__

from sim_config import load_config

CFG = load_config()


# =============================================================================
# Physical constants and defaults
#
# Values are loaded from config.yaml (see sim_config.py); this block only
# keeps the historical module-level names, since some values (R_DEFAULT) are
# derived from others and other code in this file references these names
# directly.
# =============================================================================

MU0_OVER_4PI = CFG.physics.constants.mu0_over_4pi  # T m / A
STEEL_DENSITY_DEFAULT = CFG.physics.compass_engine.steel_density  # kg/m^3
STEEL_MS_SATURATION_DEFAULT = CFG.physics.compass_engine.steel_ms_saturation  # A/m, approximately Bsat=2.0 T / mu0
# Default physical apparatus geometry used in the project report.
# The lattice parameter is defined by the centre-to-centre spacing between
# neighbouring pivots. In the code the historical variable R is kept as half
# this distance because the lattice generators use d = 2R.
CENTER_DISTANCE_DEFAULT = CFG.physics.compass_engine.center_distance  # m, centre-to-centre pivot distance
R_DEFAULT = 0.5 * CENTER_DISTANCE_DEFAULT  # m
NEEDLE_LEN_DEFAULT = CFG.physics.compass_engine.needle_len  # m, physical blade length
NEEDLE_WIDTH_DEFAULT = CFG.physics.compass_engine.needle_width  # m, physical blade width
NEEDLE_THICKNESS_DEFAULT = CFG.physics.compass_engine.needle_thickness  # m
DAMPING_DEFAULT = CFG.physics.compass_engine.damping  # N m s / rad
SOURCE_FILE_TIMESTAMP = "2026-07-10T11:50:24-03:00"  # generation/update timestamp


# =============================================================================
# Small numerical helpers
# =============================================================================


def wrap_angle(a):
    """Wrap angle(s) to (-pi, pi]. Works for plain NumPy arrays as well as
    NumPy/CuPy backend arrays (math.pi broadcasts fine against either)."""
    return (a + math.pi) % (2.0 * math.pi) - math.pi


def safe_float(x) -> float:
    """Convert NumPy/CuPy scalar to Python float."""
    try:
        if hasattr(x, "get"):
            return float(x.get())
    except Exception:
        pass
    return float(x)


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


# =============================================================================
# Physical parameter derivation
# =============================================================================


def compute_inertia_from_geometry(
    needle_len: float,
    needle_width: float,
    thickness: float,
    density: float = STEEL_DENSITY_DEFAULT,
    pivot_radius: float = 0.0,
    pivot_thickness: float = 0.0,
    pivot_density: float = 8500.0,
    pivot_mass: Optional[float] = None,
) -> float:
    """
    Moment of inertia of a rhombus-shaped needle about its central z axis.

    The blade area is approximated as 0.5 * L * W. A cylindrical pivot hole can
    be subtracted and a central pivot cylinder can be added. The pivot options
    are important for comparison with a tabletop apparatus because a metallic
    pin can contribute non-negligibly to I.
    """
    area_solid = 0.5 * needle_len * needle_width
    mass_solid = density * area_solid * thickness
    I_solid = (1.0 / 24.0) * mass_solid * (needle_len**2 + needle_width**2)

    if pivot_radius > 0.0:
        area_hole = math.pi * pivot_radius**2
        mass_hole = density * area_hole * thickness
        I_hole = 0.5 * mass_hole * pivot_radius**2
        I_blade = max(0.0, I_solid - I_hole)
    else:
        I_blade = I_solid

    if pivot_mass is not None:
        mass_pivot = pivot_mass
    else:
        mass_pivot = pivot_density * math.pi * pivot_radius**2 * pivot_thickness
    I_pivot = 0.5 * mass_pivot * pivot_radius**2
    return I_blade + I_pivot



def compute_moment_from_geometry(
    needle_len: float,
    needle_width: float,
    thickness: float,
    Ms: float = STEEL_MS_SATURATION_DEFAULT,
    pivot_radius: float = 0.0,
) -> float:
    """Magnetic moment m = Ms * volume of the rhombus blade minus pivot hole."""
    area_solid = 0.5 * needle_len * needle_width
    area_hole = math.pi * pivot_radius**2 if pivot_radius > 0.0 else 0.0
    area_net = max(0.0, area_solid - area_hole)
    return Ms * area_net * thickness


# =============================================================================
# Geometry
# =============================================================================


@dataclass
class GeometryData:
    xs: np.ndarray
    ys: np.ndarray
    theta0: np.ndarray
    r_nn: float
    Lx: float
    Ly: float
    needle_len: float
    needle_width: float
    geometry: str

    @property
    def K(self) -> int:
        return int(self.xs.size)



def make_lattice(
    N: int,
    M: int,
    geometry: str,
    R: float,
    needle_frac: float,
    noise: float,
    rng: np.random.Generator,
    needle_len: float = NEEDLE_LEN_DEFAULT,
    needle_width: float = NEEDLE_WIDTH_DEFAULT,
    use_legacy_size_from_R: bool = False,
) -> GeometryData:
    """
    Generate lattice positions and random initial angles.

    The honeycomb option follows the earlier compass.py convention: it creates
    a set of sites with hexagonal holes and returns flattened arrays shaped as
    (K,). Square and triangular arrays return K=N*M positions.

    Research default: the centre-to-centre distance is d=2R=13 mm, while the
    needle blade is an independent physical object with default dimensions
    L=10 mm and W=3 mm. This is intentionally different from the legacy
    visualization convention where needle size was set as a fraction of 2R.
    Use --use_legacy_size_from_R only for backward-compatible visual tests.
    """
    d = 2.0 * R
    s3 = math.sqrt(3.0)
    geometry = geometry.lower()

    if geometry == "square":
        jj, ii = np.meshgrid(np.arange(M), np.arange(N))
        xs = jj.astype(float) * d
        ys = ii.astype(float) * d
        Lx = M * d
        Ly = N * d

    elif geometry == "triangular":
        jj, ii = np.meshgrid(np.arange(M), np.arange(N))
        xs = jj.astype(float) * d + (ii % 2) * R
        ys = ii.astype(float) * (R * s3)
        Lx = M * d
        Ly = N * R * s3

    elif geometry == "honeycomb":
        dy = R * s3
        W = (M - 1) * d
        H = (N - 1) * dy
        n_rows = (N + 4) * 2
        x_start = -4.0 * R
        y_start = -2.0 * dy
        xs_list: List[float] = []
        ys_list: List[float] = []
        for row in range(n_rows):
            y = y_start + row * dy
            if row % 2 == 0:
                x = x_start
                while x <= W + 4.0 * R:
                    xs_list.append(x)
                    ys_list.append(y)
                    x += d
            else:
                x = x_start + R
                while x <= W + 4.0 * R:
                    xs_list.append(x)
                    ys_list.append(y)
                    x += 2.0 * d
        all_x = np.asarray(xs_list)
        all_y = np.asarray(ys_list)
        margin = 0.99 * R
        mask = (all_x >= -margin) & (all_x <= W + margin) & (all_y >= -margin) & (all_y <= H + margin)
        xs = all_x[mask]
        ys = all_y[mask]
        Lx = W
        Ly = H

    else:
        raise ValueError(f"Unknown geometry '{geometry}'. Choose square, triangular or honeycomb.")

    xs = np.asarray(xs, dtype=float).ravel()
    ys = np.asarray(ys, dtype=float).ravel()
    theta0 = noise * rng.standard_normal(xs.size)

    # Nearest-neighbour distance from the actual generated sites.
    dx = xs[:, None] - xs[None, :]
    dy = ys[:, None] - ys[None, :]
    rr = np.sqrt(dx * dx + dy * dy)
    np.fill_diagonal(rr, np.inf)
    r_nn = float(np.min(rr))

    if use_legacy_size_from_R:
        # Backward-compatible convention inherited from the visualization code:
        # changing the lattice spacing also changes the needle blade size.
        # This should not be used for calibrated research runs.
        needle_len_eff = needle_frac * d
        needle_width_eff = CFG.physics.compass_engine.legacy_width_to_length_ratio * needle_len_eff
    else:
        # Research convention: the physical needle dimensions are independent
        # of the lattice spacing. This matches the project-report apparatus.
        needle_len_eff = float(needle_len)
        needle_width_eff = float(needle_width)

    return GeometryData(xs, ys, theta0, r_nn, Lx, Ly, needle_len_eff, needle_width_eff, geometry)


# =============================================================================
# Backend and tensor precomputation
# =============================================================================


@dataclass
class Backend:
    xp: object
    name: str
    active_gpu: bool

    def asarray(self, a, dtype=None):
        return self.xp.asarray(a, dtype=dtype)

    def to_cpu(self, a):
        if self.active_gpu and hasattr(a, "get"):
            return a.get()
        return np.asarray(a)

    def sync(self):
        if self.active_gpu:
            cp.cuda.Stream.null.synchronize()  # type: ignore[union-attr]



def select_backend(use_gpu: bool) -> Backend:
    if use_gpu:
        if not _CUPY_IMPORT_OK:
            print(f"[warning] CuPy unavailable: {_CUPY_ERROR}. Falling back to CPU.")
            return Backend(np, "CPU/NumPy", False)
        try:
            cp.cuda.Device(0).compute_capability  # type: ignore[union-attr]
            test = cp.asarray([0.0, 1.0])  # type: ignore[union-attr]
            _ = cp.cos(test)  # type: ignore[union-attr]
            cp.cuda.Stream.null.synchronize()  # type: ignore[union-attr]
            return Backend(cp, "GPU/CuPy", True)  # type: ignore[arg-type]
        except Exception as exc:
            print(f"[warning] GPU initialization failed: {exc}. Falling back to CPU.")
    return Backend(np, "CPU/NumPy", False)


@dataclass
class DipolarTensor:
    Axx: object
    Axy: object
    Ayy: object
    cutoff_m: float
    cutoff_shells: float
    n_images: int
    pbc: bool



def precompute_dipolar_tensor(
    x,
    y,
    cutoff_m: float,
    cutoff_shells: float,
    backend: Backend,
    pbc: bool = False,
    Lx: Optional[float] = None,
    Ly: Optional[float] = None,
    n_images: int = 1,
    dtype=np.float64,
) -> DipolarTensor:
    """
    Precompute the linear dipolar map B_i = sum_j A_ij m_j.

    This is the main efficiency improvement. For a fixed lattice, Axx, Axy and
    Ayy do not depend on time. Time stepping then uses only matrix-vector
    products. This is much faster than rebuilding pair distances every step.
    """
    xp = backend.xp
    x = backend.asarray(x, dtype=dtype)
    y = backend.asarray(y, dtype=dtype)
    K = int(x.shape[0])

    Axx = xp.zeros((K, K), dtype=dtype)
    Axy = xp.zeros((K, K), dtype=dtype)
    Ayy = xp.zeros((K, K), dtype=dtype)

    if pbc:
        if Lx is None or Ly is None or Lx <= 0.0 or Ly <= 0.0:
            raise ValueError("PBC requires positive Lx and Ly.")
        shifts_x = [k * Lx for k in range(-n_images, n_images + 1)]
        shifts_y = [k * Ly for k in range(-n_images, n_images + 1)]
    else:
        shifts_x = [0.0]
        shifts_y = [0.0]

    for sx in shifts_x:
        for sy in shifts_y:
            rx = (x[:, None] - x[None, :]) - sx
            ry = (y[:, None] - y[None, :]) - sy
            r2 = rx * rx + ry * ry
            valid = (r2 > 1e-30) & (r2 <= cutoff_m * cutoff_m)
            r2_safe = xp.where(valid, r2, 1.0)
            r = xp.sqrt(r2_safe)
            r3 = r2_safe * r
            r5 = r2_safe * r2_safe * r
            axx = MU0_OVER_4PI * (3.0 * rx * rx / r5 - 1.0 / r3)
            axy = MU0_OVER_4PI * (3.0 * rx * ry / r5)
            ayy = MU0_OVER_4PI * (3.0 * ry * ry / r5 - 1.0 / r3)
            Axx += xp.where(valid, axx, 0.0)
            Axy += xp.where(valid, axy, 0.0)
            Ayy += xp.where(valid, ayy, 0.0)

    return DipolarTensor(Axx, Axy, Ayy, cutoff_m, cutoff_shells, n_images, pbc)


# =============================================================================
# Field protocols
# =============================================================================


@dataclass
class FieldState:
    bx: float
    by: float
    B_scalar: float
    branch: str
    curve_index: int


@dataclass
class ForcSchedule:
    starts: List[float]
    Br: List[float]
    t_down: List[float]
    t_up: List[float]
    t_sat: float
    total_time: float



def build_forc_schedule(
    Bmax: float,
    Br_min: float,
    n_curves: int,
    t_sat: float,
    t_down_default: float,
    t_up_default: float,
    rate: Optional[float],
) -> ForcSchedule:
    starts = [0.0]
    Br_list: List[float] = []
    td_list: List[float] = []
    tu_list: List[float] = []
    for k in range(n_curves):
        Br = Bmax - k * (Bmax - Br_min) / max(n_curves - 1, 1)
        span = abs(Bmax - Br)
        if rate is not None and rate > 0:
            td = span / rate
            tu = span / rate
        else:
            td = t_down_default
            tu = t_up_default
        Br_list.append(Br)
        td_list.append(td)
        tu_list.append(tu)
        starts.append(starts[-1] + t_sat + td + tu)
    return ForcSchedule(starts, Br_list, td_list, tu_list, t_sat, starts[-1])


@dataclass
class FieldProtocol:
    mode: str
    Bmax: float
    phi: float
    t_sim: float
    field_delay: float = 0.0
    t_pulse: Optional[float] = None
    sine_freq: float = 1.0
    hyst_spacing: str = "linear"
    hyst_log_k: float = 5.0
    forc: Optional[ForcSchedule] = None
    demag_freq: float = 2.0
    demag_cycles: int = 20

    @property
    def cos_phi(self) -> float:
        return math.cos(self.phi)

    @property
    def sin_phi(self) -> float:
        return math.sin(self.phi)

    def at(self, t: float) -> FieldState:
        mode = self.mode
        B = self.Bmax
        branch = mode
        idx = -1

        if mode == "static":
            Bs = B

        elif mode == "hysteresis":
            T = max(self.t_sim, 1e-30)
            t5 = T / 5.0
            if t <= t5:
                u, sgn, branch = t / t5, +1.0, "0_to_pos"
            elif t <= 2.0 * t5:
                u, sgn, branch = 1.0 - (t - t5) / t5, +1.0, "pos_to_0"
            elif t <= 3.0 * t5:
                u, sgn, branch = (t - 2.0 * t5) / t5, -1.0, "0_to_neg"
            elif t <= 4.0 * t5:
                u, sgn, branch = 1.0 - (t - 3.0 * t5) / t5, -1.0, "neg_to_0"
            else:
                u, sgn, branch = (t - 4.0 * t5) / t5, +1.0, "0_to_pos_final"
            u = max(0.0, min(1.0, u))
            if self.hyst_spacing == "log" and self.hyst_log_k > 1e-12:
                g = math.sinh(self.hyst_log_k * u) / math.sinh(self.hyst_log_k)
            else:
                g = u
            Bs = sgn * B * g

        elif mode == "forc":
            if self.forc is None:
                raise RuntimeError("FORC mode requires a schedule.")
            starts = self.forc.starts
            # Manual binary search, avoiding scipy dependency.
            k = np.searchsorted(starts, t, side="right") - 1
            k = int(max(0, min(k, len(self.forc.Br) - 1)))
            idx = k
            tl = t - starts[k]
            ts = self.forc.t_sat
            td = self.forc.t_down[k]
            tu = self.forc.t_up[k]
            Br = self.forc.Br[k]
            if tl <= ts:
                Bs = B
                branch = "sat"
            elif tl <= ts + td:
                f = (tl - ts) / max(td, 1e-30)
                Bs = (1.0 - f) * B + f * Br
                branch = "down"
            else:
                f = (tl - ts - td) / max(tu, 1e-30)
                Bs = (1.0 - f) * Br + f * B
                branch = "up"

        elif mode == "sine":
            Bs = B * math.sin(2.0 * math.pi * self.sine_freq * t)
            branch = "sine"

        elif mode in ("step_up", "step_pos"):
            Bs = 0.0 if t < self.field_delay else B
            branch = "before_step" if t < self.field_delay else "after_step"

        elif mode in ("step_down", "step_neg"):
            Bs = B if t < self.field_delay else 0.0
            branch = "before_step" if t < self.field_delay else "after_step"

        elif mode == "pulse":
            if t < self.field_delay:
                Bs = 0.0
                branch = "delay"
            elif self.t_pulse is not None and t >= self.field_delay + self.t_pulse:
                Bs = 0.0
                branch = "post_pulse"
            else:
                Bs = B
                branch = "pulse_on"

        elif mode == "demag_rot":
            Tdemag = self.demag_cycles / max(self.demag_freq, 1e-30)
            if t <= Tdemag:
                envelope = B * (1.0 - t / max(Tdemag, 1e-30))
                angle = 2.0 * math.pi * self.demag_freq * t
                return FieldState(envelope * math.cos(angle), envelope * math.sin(angle), envelope, "demag_rot", -1)
            Bs = 0.0
            branch = "post_demag"

        elif mode == "demag_linear":
            Tdemag = self.demag_cycles / max(self.demag_freq, 1e-30)
            Bs = B * max(0.0, 1.0 - t / max(Tdemag, 1e-30))
            branch = "demag_linear" if t <= Tdemag else "post_demag"

        else:
            raise ValueError(f"Unknown field_mode '{mode}'.")

        return FieldState(Bs * self.cos_phi, Bs * self.sin_phi, Bs, branch, idx)


# =============================================================================
# Metrics and domain statistics
# =============================================================================


@dataclass
class Metrics:
    Mx: float
    My: float
    M_proj: float
    S1: float
    S2: float
    theta_director: float
    q_axis: float
    omega_rms: float
    omega_max: float
    E_dip: float
    E_ext: float
    E_kin: float



def compute_metrics(
    theta,
    omega,
    tensor: DipolarTensor,
    moment: float,
    inertia: float,
    bx: float,
    by: float,
    phi_drive: float,
    backend: Backend,
) -> Metrics:
    xp = backend.xp
    mx = moment * xp.cos(theta)
    my = moment * xp.sin(theta)

    Bdx = tensor.Axx @ mx + tensor.Axy @ my
    Bdy = tensor.Axy @ mx + tensor.Ayy @ my

    Mx = xp.mean(xp.cos(theta))
    My = xp.mean(xp.sin(theta))
    M_proj = Mx * math.cos(phi_drive) + My * math.sin(phi_drive)
    S1 = xp.sqrt(Mx * Mx + My * My)

    c2 = xp.mean(xp.cos(2.0 * theta))
    s2 = xp.mean(xp.sin(2.0 * theta))
    S2 = xp.sqrt(c2 * c2 + s2 * s2)
    theta_director = 0.5 * xp.arctan2(s2, c2)

    # q_axis measures population imbalance between x-like and y-like axes.
    # It is meaningful for square-type director studies, but still useful as a
    # diagnostic in other geometries if interpreted cautiously.
    x_like = xp.abs(xp.cos(theta)) >= xp.abs(xp.sin(theta))
    q_axis = xp.abs(xp.mean(xp.where(x_like, 1.0, -1.0)))

    E_dip = -0.5 * xp.sum(mx * Bdx + my * Bdy)
    E_ext = -xp.sum(mx * bx + my * by)
    E_kin = 0.5 * inertia * xp.sum(omega * omega)
    omega_rms = xp.sqrt(xp.mean(omega * omega))
    omega_max = xp.max(xp.abs(omega))

    return Metrics(
        safe_float(Mx),
        safe_float(My),
        safe_float(M_proj),
        safe_float(S1),
        safe_float(S2),
        safe_float(theta_director),
        safe_float(q_axis),
        safe_float(omega_rms),
        safe_float(omega_max),
        safe_float(E_dip),
        safe_float(E_ext),
        safe_float(E_kin),
    )



def domain_statistics(theta_cpu: np.ndarray, x_cpu: np.ndarray, y_cpu: np.ndarray, r_nn: float, tol_deg: float) -> Dict[str, object]:
    """
    Simple connected-component domain statistics.

    Two needles are in the same local domain when they are spatial nearest
    neighbours and their angles differ by less than tol_deg modulo 2pi. This is
    a strict polar-domain definition. For director domains, users should repeat
    the analysis with theta modulo pi in a separate script.
    """
    K = int(theta_cpu.size)
    if K == 0:
        return {"n_domains": 0, "mean_domain_size": 0.0, "max_domain_size": 0, "domain_sizes": []}
    parent = np.arange(K)
    size = np.ones(K, dtype=int)

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a: int, b: int):
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if size[ra] < size[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        size[ra] += size[rb]

    tol = math.radians(tol_deg)
    # O(K^2), but called only at the end. Avoids scipy dependency.
    for i in range(K):
        dx = x_cpu[i] - x_cpu[i + 1 :]
        dy = y_cpu[i] - y_cpu[i + 1 :]
        dist = np.sqrt(dx * dx + dy * dy)
        neigh = np.where(dist <= 1.05 * r_nn)[0] + i + 1
        if neigh.size == 0:
            continue
        dtheta = np.abs(wrap_angle(theta_cpu[i] - theta_cpu[neigh]))
        for j in neigh[dtheta <= tol]:
            union(i, int(j))

    roots = np.array([find(i) for i in range(K)])
    _, counts = np.unique(roots, return_counts=True)
    counts = np.sort(counts)[::-1]
    return {
        "n_domains": int(counts.size),
        "mean_domain_size": float(np.mean(counts)),
        "max_domain_size": int(counts[0]),
        "domain_sizes": [int(c) for c in counts.tolist()],
    }


# =============================================================================
# Simulation engine
# =============================================================================


@dataclass
class SimulationConfig:
    geometry: str
    N: int
    M: int
    R: float
    center_distance: float
    needle_frac: float
    needle_len: float
    needle_width: float
    use_legacy_size_from_R: bool
    needle_thickness: float
    steel_density: float
    steel_Ms: float
    moment: float
    inertia: float
    damping: float
    damping_noise: float
    t_sim: float
    dt_factor: float
    field_mode: str
    B_ext: float
    phi_ext_deg: float
    cutoff_m: float
    cutoff_shells: float
    pbc: bool
    n_images: int
    seed: int
    backend: str
    log_every: int
    flip_angle_deg: float
    flip_band_deg: float
    flip_dwell_T0: float
    flip_settle_frac: float
    dt_guard_alpha: float
    dt_guard_substep: bool
    event_log: bool
    domain_tol_deg: float



def run_simulation(args) -> Tuple[Path, Path, Path]:
    out_dir = ensure_dir(args.out_dir)
    data_dir = ensure_dir(out_dir / "data")
    meta_dir = ensure_dir(out_dir / "meta")
    state_dir = ensure_dir(out_dir / "states")
    image_dir = ensure_dir(out_dir / "images")

    seed = args.seed if args.seed is not None else int(time.time_ns() % (2**32 - 1))
    rng = np.random.default_rng(seed)
    geom = make_lattice(
        args.N,
        args.M,
        args.geometry,
        args.R,
        args.needle_frac,
        args.noise,
        rng,
        needle_len=args.needle_len,
        needle_width=args.needle_width,
        use_legacy_size_from_R=bool(args.use_legacy_size_from_R),
    )

    needle_thickness = args.needle_thickness
    steel_Ms = args.steel_Ms
    if args.steel_Bsat is not None:
        # Bsat = mu0 Ms, so Ms = Bsat / mu0.
        steel_Ms = args.steel_Bsat / (4.0 * math.pi * 1.0e-7)

    moment = args.moment if args.moment is not None else compute_moment_from_geometry(
        geom.needle_len,
        geom.needle_width,
        needle_thickness,
        steel_Ms,
        args.pivot_radius,
    )
    inertia = args.inertia if args.inertia is not None else compute_inertia_from_geometry(
        geom.needle_len,
        geom.needle_width,
        needle_thickness,
        args.steel_density,
        args.pivot_radius,
        args.pivot_thickness,
        args.pivot_density,
        args.pivot_mass,
    )

    B_ref = MU0_OVER_4PI * 2.0 * moment / geom.r_nn**3
    B_ext = args.B_ext if args.B_ext is not None else args.B_max_factor * B_ref
    phi = math.radians(args.phi_ext_deg)

    if args.cutoff_m is not None:
        cutoff_m = float(args.cutoff_m)
        cutoff_shells = cutoff_m / geom.r_nn
        cutoff_policy = "absolute_m"
    else:
        cutoff_shells = float(args.cutoff_shells)
        cutoff_m = cutoff_shells * geom.r_nn
        cutoff_policy = "shells_times_r_nn"

    B_eff = max(abs(B_ext), B_ref)
    omega0 = math.sqrt(moment * B_eff / inertia)
    T0 = 2.0 * math.pi / omega0
    dt = args.dt_factor * T0

    if args.field_mode == "forc":
        forc_Br_min = args.forc_Br_min if args.forc_Br_min is not None else -abs(B_ext)
        forc_schedule = build_forc_schedule(
            abs(B_ext),
            forc_Br_min,
            args.forc_n_curves,
            args.forc_t_sat,
            args.forc_t_ramp_down,
            args.forc_t_ramp_up,
            args.forc_rate,
        )
        t_sim = forc_schedule.total_time
    elif args.field_mode == "demag_rot" or args.field_mode == "demag_linear":
        forc_schedule = None
        t_sim = args.demag_cycles / max(args.demag_freq, 1e-30) + args.t_relax_after
    else:
        forc_schedule = None
        t_sim = args.t_sim

    protocol = FieldProtocol(
        mode=args.field_mode,
        Bmax=B_ext,
        phi=phi,
        t_sim=t_sim,
        field_delay=args.field_delay,
        t_pulse=args.t_pulse,
        sine_freq=args.field_freq,
        hyst_spacing=args.hyst_spacing,
        hyst_log_k=args.hyst_log_k,
        forc=forc_schedule,
        demag_freq=args.demag_freq,
        demag_cycles=args.demag_cycles,
    )

    backend = select_backend(args.use_gpu)
    dtype = np.float32 if args.float32 else np.float64
    tensor_bytes = 3.0 * geom.K * geom.K * np.dtype(dtype).itemsize
    if tensor_bytes > args.tensor_mem_limit_gb * 1e9:
        raise MemoryError(
            f"Dipolar tensor would require {tensor_bytes/1e9:.2f} GB. "
            f"Increase --tensor_mem_limit_gb or reduce grid size."
        )

    tensor = precompute_dipolar_tensor(
        geom.xs,
        geom.ys,
        cutoff_m,
        cutoff_shells,
        backend,
        pbc=args.pbc,
        Lx=geom.Lx,
        Ly=geom.Ly,
        n_images=args.n_images,
        dtype=dtype,
    )

    xp = backend.xp
    theta = backend.asarray(geom.theta0, dtype=dtype)
    omega = xp.zeros(geom.K, dtype=dtype)

    if args.damping_noise > 0.0:
        damping_np = args.damping * (1.0 + args.damping_noise * rng.uniform(-1.0, 1.0, geom.K))
        damping_np = np.maximum(damping_np, 0.0)
        damping = backend.asarray(damping_np, dtype=dtype)
    else:
        damping = float(args.damping)

    Q = omega0 * inertia / args.damping if args.damping > 0 else math.inf
    n_steps = max(1, int(math.ceil(t_sim / dt)))
    log_every = max(1, args.log_every)
    flip_threshold = math.radians(args.flip_angle_deg)

    # ------------------------------------------------------------------
    # V2.1 hardened flip detection: derived quantities.
    #
    # Drive-axis channel (flip_field): a needle is in state +1 when its
    # projection cos(theta - phi) >= +u_thresh, in state -1 when it is
    # <= -u_thresh, and "undecided" inside the dead band. The dead band is
    # angular: half-width flip_band_deg measured FROM THE PERPENDICULAR to
    # the drive axis, so u_thresh = sin(flip_band_deg). A state change is
    # committed only after the candidate state has persisted for
    # n_dwell_steps consecutive steps.
    #
    # Rest-angle channel (flip_angle): each needle carries a reference
    # angle theta_ref (its last settled orientation). An event is
    # committed when |wrap(theta - theta_ref)| >= flip_threshold AND
    # |omega| <= omega_settle, sustained for n_dwell_steps. theta_ref is
    # then updated. This channel needs no knowledge of the local energy
    # landscape and works identically for all lattice geometries.
    # ------------------------------------------------------------------
    u_thresh = math.sin(math.radians(args.flip_band_deg))
    n_dwell_steps = max(1, int(round(args.flip_dwell_T0 * T0 / dt)))
    omega_settle = args.flip_settle_frac * omega0

    # Initial field and torque.
    f0 = protocol.at(0.0)
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    def torque_and_field(th, bx, by):
        """Return (torque, Bx_total, By_total) per needle.

        The total field components are returned because they are already
        computed by the matrix-vector products; the stability monitor uses
        them at zero additional matvec cost.
        """
        mx = moment * xp.cos(th)
        my = moment * xp.sin(th)
        Bdx = tensor.Axx @ mx + tensor.Axy @ my
        Bdy = tensor.Axy @ mx + tensor.Ayy @ my
        Btx = Bdx + bx
        Bty = Bdy + by
        return mx * Bty - my * Btx, Btx, Bty

    tau, Btx, Bty = torque_and_field(theta, f0.bx, f0.by)

    # Drive-axis Schmitt state. Needles initially inside the dead band get
    # state 0 (undecided); their first commitment is not counted as a flip.
    u0 = xp.cos(theta - phi)
    sigma_state = xp.where(u0 >= u_thresh, 1, xp.where(u0 <= -u_thresh, -1, 0)).astype(xp.int8)
    ff_pending = xp.zeros(geom.K, dtype=xp.int8)
    ff_count = xp.zeros(geom.K, dtype=xp.int32)

    # Rest-angle channel state.
    theta_ref = theta.copy()
    fa_count = xp.zeros(geom.K, dtype=xp.int32)

    flip_field_acc = 0
    flip_angle_acc = 0
    flip_field_total = 0
    flip_angle_total = 0

    # Event log buffer: (step, t, needle_id, channel, theta_committed).
    # event_path is assigned after the run tag is constructed below.
    event_buffer: List[Tuple[int, float, int, str, float]] = []
    event_path: Optional[Path] = None

    # Stability monitor.
    guard_alpha = float(args.dt_guard_alpha)
    guard_flagged_steps = 0
    guard_substepped_steps = 0
    guard_max_ratio = 0.0

    invI = 1.0 / inertia

    tag = args.tag or f"{args.geometry}_{args.field_mode}_N{args.N}_M{args.M}_seed{seed}"
    if args.event_log:
        event_path = data_dir / f"{tag}_events.csv"
    csv_path = data_dir / f"{tag}.csv"
    meta_path = meta_dir / f"{tag}.json"
    initial_state_path = state_dir / f"{tag}_initial.npz"
    state_path = state_dir / f"{tag}_final.npz"

    config = SimulationConfig(
        geometry=args.geometry,
        N=args.N,
        M=args.M,
        R=args.R,
        center_distance=2.0 * args.R,
        needle_frac=args.needle_frac,
        needle_len=geom.needle_len,
        needle_width=geom.needle_width,
        use_legacy_size_from_R=bool(args.use_legacy_size_from_R),
        needle_thickness=needle_thickness,
        steel_density=args.steel_density,
        steel_Ms=steel_Ms,
        moment=moment,
        inertia=inertia,
        damping=args.damping,
        damping_noise=args.damping_noise,
        t_sim=t_sim,
        dt_factor=args.dt_factor,
        field_mode=args.field_mode,
        B_ext=B_ext,
        phi_ext_deg=args.phi_ext_deg,
        cutoff_m=cutoff_m,
        cutoff_shells=cutoff_shells,
        pbc=args.pbc,
        n_images=args.n_images,
        seed=seed,
        backend=backend.name,
        log_every=log_every,
        flip_angle_deg=args.flip_angle_deg,
        flip_band_deg=args.flip_band_deg,
        flip_dwell_T0=args.flip_dwell_T0,
        flip_settle_frac=args.flip_settle_frac,
        dt_guard_alpha=args.dt_guard_alpha,
        dt_guard_substep=bool(args.dt_guard_substep),
        event_log=bool(args.event_log),
        domain_tol_deg=args.domain_tol_deg,
    )

    # Nominal drive sweep rate, logged so rate can be a first-class campaign
    # axis (athermal deterministic dynamics is rate-dependent by construction).
    if args.field_mode == "hysteresis":
        sweep_rate_T_s = abs(B_ext) / max(t_sim / 5.0, 1e-30)
    elif args.field_mode == "forc" and args.forc_rate is not None:
        sweep_rate_T_s = float(args.forc_rate)
    else:
        sweep_rate_T_s = None

    metadata = {
        "program": "compassV02.py",
        "version": "2.1.0",
        "source_file_timestamp": SOURCE_FILE_TIMESTAMP,
        "created_unix_time": time.time(),
        "created_datetime_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "created_datetime_local": datetime.now().astimezone().isoformat(timespec="seconds"),
        "tag": tag,
        "config": asdict(config),
        "derived": {
            "K": geom.K,
            "r_nn_m": geom.r_nn,
            "center_distance_default_m": CENTER_DISTANCE_DEFAULT,
            "Lx_m": geom.Lx,
            "Ly_m": geom.Ly,
            "needle_len_m": geom.needle_len,
            "needle_width_m": geom.needle_width,
            "B_ref_T": B_ref,
            "B_eff_T": B_eff,
            "omega0_rad_s": omega0,
            "T0_s": T0,
            "dt_s": dt,
            "n_steps": n_steps,
            "Q": Q,
            "cutoff_policy": cutoff_policy,
            "tensor_bytes": tensor_bytes,
            "dtype": "float32" if args.float32 else "float64",
            "u_thresh": u_thresh,
            "n_dwell_steps": n_dwell_steps,
            "omega_settle_rad_s": omega_settle,
            "sweep_rate_T_per_s": sweep_rate_T_s,
            "sweep_rate_Bref_per_T0": (sweep_rate_T_s * T0 / B_ref) if (sweep_rate_T_s and B_ref > 0.0) else None,
        },
        "notes": {
            "S1": "polar order |<exp(i theta)>|",
            "S2": "nematic/director order |<exp(2 i theta)>|",
            "flip_field": "COMMITTED drive-axis reversals since previous log row (Schmitt band flip_band_deg + dwell flip_dwell_T0); ringing-proof",
            "flip_angle": "COMMITTED rest-angle displacement events since previous log row (>= flip_angle_deg from last settled angle, |omega| settled); geometry-agnostic",
            "E_dip": "-0.5 sum_i m_i dot B_i^dip",
            "E_ext": "-sum_i m_i dot B_ext",
            "E_kin": "0.5 I sum_i omega_i^2",
        },
    }

    # Save the initial state before time integration.
    np.savez_compressed(
        initial_state_path,
        xs=geom.xs,
        ys=geom.ys,
        theta=backend.to_cpu(theta),
        omega=backend.to_cpu(omega),
        r_nn=geom.r_nn,
        metadata_json=json.dumps(metadata),
    )

    if args.verbose:
        print(f"Backend      : {backend.name}")
        print(f"Geometry     : {args.geometry}, K={geom.K}, r_nn={geom.r_nn:.6g} m")
        print(f"Needle       : L={geom.needle_len*1e3:.3f} mm, W={geom.needle_width*1e3:.3f} mm, t={needle_thickness*1e3:.3f} mm")
        print(f"Center dist. : {2.0*args.R*1e3:.3f} mm (R={args.R*1e3:.3f} mm)")
        print(f"B_ref        : {B_ref*1e3:.6g} mT")
        print(f"B_ext        : {B_ext*1e3:.6g} mT")
        print(f"omega0, T0   : {omega0:.6g} rad/s, {T0:.6g} s")
        print(f"dt, n_steps  : {dt:.6g} s, {n_steps}")
        print(f"Q            : {Q:.6g}")
        print(f"cutoff       : {cutoff_m:.6g} m = {cutoff_shells:.3g} r_nn")
        print(f"tensor       : {tensor_bytes/1e6:.1f} MB")
        print(f"output       : {csv_path}")

    fieldnames = [
        "step",
        "t_s",
        "Bx_T",
        "By_T",
        "B_scalar_T",
        "branch",
        "forc_index",
        "Mx",
        "My",
        "M_proj",
        "S1",
        "S2",
        "theta_director_rad",
        "q_axis",
        "flip_field",
        "flip_angle",
        "E_dip_J",
        "E_ext_J",
        "E_kin_J",
        "E_total_J",
        "omega_rms_rad_s",
        "omega_max_rad_s",
    ]

    t_start = time.perf_counter()
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        def write_row(step: int, t: float, f: FieldState):
            nonlocal flip_field_acc, flip_angle_acc
            met = compute_metrics(theta, omega, tensor, moment, inertia, f.bx, f.by, phi, backend)
            writer.writerow(
                {
                    "step": step,
                    "t_s": f"{t:.12g}",
                    "Bx_T": f"{f.bx:.12g}",
                    "By_T": f"{f.by:.12g}",
                    "B_scalar_T": f"{f.B_scalar:.12g}",
                    "branch": f.branch,
                    "forc_index": f.curve_index,
                    "Mx": f"{met.Mx:.12g}",
                    "My": f"{met.My:.12g}",
                    "M_proj": f"{met.M_proj:.12g}",
                    "S1": f"{met.S1:.12g}",
                    "S2": f"{met.S2:.12g}",
                    "theta_director_rad": f"{met.theta_director:.12g}",
                    "q_axis": f"{met.q_axis:.12g}",
                    "flip_field": int(flip_field_acc),
                    "flip_angle": int(flip_angle_acc),
                    "E_dip_J": f"{met.E_dip:.12g}",
                    "E_ext_J": f"{met.E_ext:.12g}",
                    "E_kin_J": f"{met.E_kin:.12g}",
                    "E_total_J": f"{(met.E_dip + met.E_ext + met.E_kin):.12g}",
                    "omega_rms_rad_s": f"{met.omega_rms:.12g}",
                    "omega_max_rad_s": f"{met.omega_max:.12g}",
                }
            )
            flip_field_acc = 0
            flip_angle_acc = 0

        write_row(0, 0.0, f0)

        def verlet_step(th, om, ta, bx, by, dt_s):
            """One damped velocity-Verlet step of size dt_s. Returns new state."""
            bh = damping * dt_s * 0.5 * invI
            c_om = (1.0 - bh) / (1.0 + bh)
            c_ta = dt_s * 0.5 * invI / (1.0 + bh)
            acc = (ta - damping * om) * invI
            th_n = th + om * dt_s + 0.5 * dt_s * dt_s * acc
            th_n = wrap_angle(th_n)
            ta_n, btx, bty = torque_and_field(th_n, bx, by)
            om_n = om * c_om + (ta + ta_n) * c_ta
            return th_n, om_n, ta_n, btx, bty

        def flush_events():
            if event_path is None or not event_buffer:
                return
            new_file = not event_path.exists()
            with open(event_path, "a", newline="") as efh:
                ew = csv.writer(efh)
                if new_file:
                    ew.writerow(["step", "t_s", "needle_id", "channel", "theta_rad"])
                ew.writerows(event_buffer)
            event_buffer.clear()

        for step in range(1, n_steps + 1):
            t = step * dt
            f = protocol.at(t)

            theta_new, omega_new, tau_new, Btx, Bty = verlet_step(theta, omega, tau, f.bx, f.by, dt)

            # --------------------------------------------------------------
            # Stability monitor: local stiffness frequency from |B_i|.
            # tau_i = (m x B_i)_z, so the linearized local frequency is
            # bounded by sqrt(m |B_i| / I). Uses fields already computed.
            # --------------------------------------------------------------
            B2max = safe_float(xp.max(Btx * Btx + Bty * Bty))
            ratio = math.sqrt(moment * math.sqrt(B2max) * invI) * dt
            if ratio > guard_max_ratio:
                guard_max_ratio = ratio
            if ratio > guard_alpha:
                guard_flagged_steps += 1
                if args.dt_guard_substep:
                    # Re-integrate this step globally with 4 sub-steps.
                    # A partial (per-needle) substep is not an option: the
                    # dipolar field is global and would desynchronize.
                    guard_substepped_steps += 1
                    th_s, om_s, ta_s = theta, omega, tau
                    dt_s = 0.25 * dt
                    for k_sub in range(4):
                        f_s = protocol.at((step - 1) * dt + (k_sub + 1) * dt_s)
                        th_s, om_s, ta_s, Btx, Bty = verlet_step(th_s, om_s, ta_s, f_s.bx, f_s.by, dt_s)
                    theta_new, omega_new, tau_new = th_s, om_s, ta_s

            # --------------------------------------------------------------
            # Channel A: drive-axis reversal (flip_field), Schmitt + dwell.
            # --------------------------------------------------------------
            u = xp.cos(theta_new - phi)
            cand = xp.where(u >= u_thresh, 1, xp.where(u <= -u_thresh, -1, 0)).astype(xp.int8)
            in_new_state = (cand != 0) & (cand != sigma_state)
            same_as_pending = in_new_state & (cand == ff_pending)
            ff_count = xp.where(same_as_pending, ff_count + 1,
                                xp.where(in_new_state, 1, 0)).astype(xp.int32)
            ff_pending = xp.where(in_new_state, cand, 0).astype(xp.int8)
            commit_ff = ff_count >= n_dwell_steps
            n_commit_ff = int(safe_float(xp.sum(commit_ff)))
            if n_commit_ff > 0:
                # First commitment of initially undecided needles (state 0)
                # establishes the state without counting as a reversal.
                real_flip = commit_ff & (sigma_state != 0)
                n_real = int(safe_float(xp.sum(real_flip)))
                flip_field_acc += n_real
                flip_field_total += n_real
                if event_path is not None and n_real > 0:
                    ids = backend.to_cpu(xp.where(real_flip)[0])
                    ths = backend.to_cpu(theta_new)
                    event_buffer.extend(
                        (step, t, int(i), "field", float(ths[int(i)])) for i in ids
                    )
                sigma_state = xp.where(commit_ff, ff_pending, sigma_state).astype(xp.int8)
                ff_count = xp.where(commit_ff, 0, ff_count).astype(xp.int32)
                ff_pending = xp.where(commit_ff, 0, ff_pending).astype(xp.int8)

            # --------------------------------------------------------------
            # Channel B: rest-angle displacement (flip_angle), settle + dwell.
            # --------------------------------------------------------------
            disp = xp.abs(wrap_angle(theta_new - theta_ref)) >= flip_threshold
            settled = xp.abs(omega_new) <= omega_settle
            cond = disp & settled
            fa_count = xp.where(cond, fa_count + 1, 0).astype(xp.int32)
            commit_fa = fa_count >= n_dwell_steps
            n_commit_fa = int(safe_float(xp.sum(commit_fa)))
            if n_commit_fa > 0:
                flip_angle_acc += n_commit_fa
                flip_angle_total += n_commit_fa
                if event_path is not None:
                    ids = backend.to_cpu(xp.where(commit_fa)[0])
                    ths = backend.to_cpu(theta_new)
                    event_buffer.extend(
                        (step, t, int(i), "angle", float(ths[int(i)])) for i in ids
                    )
                theta_ref = xp.where(commit_fa, theta_new, theta_ref)
                fa_count = xp.where(commit_fa, 0, fa_count).astype(xp.int32)

            theta = theta_new
            omega = omega_new
            tau = tau_new

            if step % log_every == 0 or step == n_steps:
                write_row(step, t, f)
                if len(event_buffer) >= 4096:
                    flush_events()

            if args.progress and step % max(1, n_steps // 100) == 0:
                elapsed = time.perf_counter() - t_start
                rate = step / max(elapsed, 1e-12)
                print(f"\r{100*step/n_steps:6.2f}%  step {step}/{n_steps}  {rate:8.0f} steps/s", end="", flush=True)

        flush_events()

    if args.progress:
        print()
    backend.sync()

    theta_cpu = backend.to_cpu(theta)
    omega_cpu = backend.to_cpu(omega)
    domain_stats = domain_statistics(theta_cpu, geom.xs, geom.ys, geom.r_nn, args.domain_tol_deg)
    metadata["domain_statistics_final"] = domain_stats
    metadata["runtime_s"] = time.perf_counter() - t_start
    metadata["flip_statistics"] = {
        "flip_field_total_committed": int(flip_field_total),
        "flip_angle_total_committed": int(flip_angle_total),
        "u_thresh": u_thresh,
        "n_dwell_steps": int(n_dwell_steps),
        "omega_settle_rad_s": omega_settle,
    }
    metadata["stability_monitor"] = {
        "dt_guard_alpha": guard_alpha,
        "flagged_steps": int(guard_flagged_steps),
        "substepped_steps": int(guard_substepped_steps),
        "max_omega_local_dt": guard_max_ratio,
        "note": "omega_local = sqrt(m|B_i|/I); flagged when max_i(omega_local)*dt > alpha",
    }
    if event_path is not None:
        metadata["event_log_path"] = str(event_path)
    if guard_flagged_steps > 0 and not args.dt_guard_substep:
        print(f"[warning] stability monitor flagged {guard_flagged_steps} steps "
              f"(max omega_local*dt = {guard_max_ratio:.3f} > alpha = {guard_alpha}). "
              f"Consider smaller --dt_factor or --dt_guard_substep, and record in quality control.")

    np.savez_compressed(
        state_path,
        xs=geom.xs,
        ys=geom.ys,
        theta=theta_cpu,
        omega=omega_cpu,
        r_nn=geom.r_nn,
        metadata_json=json.dumps(metadata),
    )

    with open(meta_path, "w") as fh:
        json.dump(metadata, fh, indent=2)

    # Always generate initial and final lattice figures.
    try:
        make_lattice_png(
            initial_state_path,
            image_dir / f"{tag}_initial_lattice.png",
            dpi=args.png_dpi,
            transparent=args.png_transparent,
            with_axes=args.png_with_axes,
            title="Initial state",
        )
        make_lattice_png(
            state_path,
            image_dir / f"{tag}_final_lattice.png",
            dpi=args.png_dpi,
            transparent=args.png_transparent,
            with_axes=args.png_with_axes,
            title="Final state",
        )
        make_initial_final_lattice_png(
            initial_state_path,
            state_path,
            image_dir / f"{tag}_initial_final.png",
            dpi=args.png_dpi,
            transparent=args.png_transparent,
            with_axes=args.png_with_axes,
            panel_titles=not args.png_no_panel_titles,
        )
    except Exception as exc:
        print(f"[warning] lattice PNG generation failed: {exc}")

    if args.make_plot:
        try:
            make_quick_plot(csv_path, out_dir / f"{tag}_quicklook.png")
        except Exception as exc:
            print(f"[warning] quick plot failed: {exc}")

    return csv_path, meta_path, state_path



# =============================================================================
# Lattice PNG plotting
# =============================================================================


def _load_state_npz_for_png(npz_path: Path) -> Dict[str, object]:
    data = np.load(npz_path, allow_pickle=True)
    meta = json.loads(str(data["metadata_json"])) if "metadata_json" in data.files else {}
    derived = meta.get("derived", {})
    return {
        "xs": np.asarray(data["xs"], dtype=float),
        "ys": np.asarray(data["ys"], dtype=float),
        "theta": np.asarray(data["theta"], dtype=float),
        "r_nn": float(data["r_nn"]) if "r_nn" in data.files else float(derived.get("r_nn_m", 1.0)),
        "needle_len": float(derived.get("needle_len_m", NEEDLE_LEN_DEFAULT)),
        "needle_width": float(derived.get("needle_width_m", NEEDLE_WIDTH_DEFAULT)),
        "metadata": meta,
    }


def _needle_halves_for_png(x: float, y: float, theta: float, length: float, width: float) -> Tuple[np.ndarray, np.ndarray]:
    """Return blue/north and white/south triangular halves of a compass needle."""
    u = np.array([math.cos(theta), math.sin(theta)])
    v = np.array([-math.sin(theta), math.cos(theta)])
    c = np.array([x, y])
    left_tip = c - 0.5 * length * u
    right_tip = c + 0.5 * length * u
    top_mid = c + 0.5 * width * v
    bottom_mid = c - 0.5 * width * v
    return np.array([left_tip, top_mid, bottom_mid]), np.array([right_tip, top_mid, bottom_mid])


def _draw_lattice_png(ax, state: Dict[str, object], with_axes: bool = False, title: Optional[str] = None):
    from matplotlib.patches import Polygon, Circle

    xs = state["xs"]
    ys = state["ys"]
    theta = state["theta"]
    r_nn = float(state["r_nn"])
    needle_len = float(state["needle_len"])
    needle_width = float(state["needle_width"])

    colors = CFG.physics.needle_render.colors
    blue_north = colors["blue_north"]
    white_south = colors["white_south"]
    edge = colors["edge"]
    pivot = colors["pivot"]
    pivot_radius_frac = CFG.physics.needle_render.pivot_radius_frac
    pivot_inner_radius_frac = CFG.physics.needle_render.pivot_inner_radius_frac

    for x, y, th in zip(xs, ys, theta):
        north_half, south_half = _needle_halves_for_png(float(x), float(y), float(th), needle_len, needle_width)
        ax.add_patch(Polygon(north_half, closed=True, facecolor=blue_north, edgecolor=edge,
                             linewidth=0.8, joinstyle="miter", zorder=2))
        ax.add_patch(Polygon(south_half, closed=True, facecolor=white_south, edgecolor=edge,
                             linewidth=0.8, joinstyle="miter", zorder=2))
        ax.add_patch(Circle((float(x), float(y)), pivot_radius_frac * r_nn, facecolor=pivot,
                            edgecolor=edge, linewidth=0.6, zorder=5))
        ax.add_patch(Circle((float(x), float(y)), pivot_inner_radius_frac * r_nn, facecolor=pivot,
                            edgecolor="white", linewidth=0.35, zorder=6))

    margin = 0.8 * r_nn
    ax.set_xlim(float(np.min(xs)) - margin, float(np.max(xs)) + margin)
    ax.set_ylim(float(np.min(ys)) - margin, float(np.max(ys)) + margin)
    ax.set_aspect("equal")

    if with_axes:
        if title:
            ax.set_title(title)
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.grid(alpha=0.25)
    else:
        ax.axis("off")


def make_lattice_png(npz_path: Path, png_path: Path, dpi: int = 300,
                     transparent: bool = False, with_axes: bool = False,
                     title: Optional[str] = None):
    """Generate one clean lattice PNG from a saved state NPZ."""
    import matplotlib.pyplot as plt

    state = _load_state_npz_for_png(npz_path)
    fig, ax = plt.subplots(figsize=(6, 6), facecolor="none" if transparent else "white")
    _draw_lattice_png(ax, state, with_axes=with_axes, title=title)
    if with_axes:
        fig.tight_layout()
    else:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ensure_dir(Path(png_path).parent)
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight", pad_inches=0 if not with_axes else 0.1,
                facecolor="none" if transparent else "white", transparent=transparent)
    plt.close(fig)


def make_initial_final_lattice_png(initial_npz: Path, final_npz: Path, png_path: Path,
                                   dpi: int = 300, transparent: bool = False,
                                   with_axes: bool = False, panel_titles: bool = True):
    """Generate side-by-side initial/final lattice PNG from saved state NPZ files."""
    import matplotlib.pyplot as plt

    s0 = _load_state_npz_for_png(initial_npz)
    s1 = _load_state_npz_for_png(final_npz)
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), facecolor="none" if transparent else "white")
    _draw_lattice_png(axes[0], s0, with_axes=with_axes, title="Initial state" if panel_titles else None)
    _draw_lattice_png(axes[1], s1, with_axes=with_axes, title="Final state" if panel_titles else None)
    if with_axes:
        fig.tight_layout()
    else:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0.02)
    ensure_dir(Path(png_path).parent)
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight", pad_inches=0 if not with_axes else 0.1,
                facecolor="none" if transparent else "white", transparent=transparent)
    plt.close(fig)

# =============================================================================
# Quick plotting
# =============================================================================



def make_quick_plot(csv_path: Path, png_path: Path):
    """Simple diagnostic plot. Kept separate so the simulation engine has no plotting dependency at runtime unless requested."""
    import matplotlib.pyplot as plt

    data = np.genfromtxt(csv_path, delimiter=",", names=True, dtype=None, encoding=None)
    B = np.asarray(data["B_scalar_T"], dtype=float) * 1e3
    M = np.asarray(data["M_proj"], dtype=float)
    t = np.asarray(data["t_s"], dtype=float)
    S1 = np.asarray(data["S1"], dtype=float)
    S2 = np.asarray(data["S2"], dtype=float)

    fig, ax = plt.subplots(figsize=(6.5, 5.0))
    ax.plot(B, M, lw=1.2)
    ax.set_xlabel("B along drive axis (mT)")
    ax.set_ylabel("M projection")
    ax.set_title("M(B)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(png_path, dpi=160)
    plt.close(fig)

    fig2, ax2 = plt.subplots(figsize=(7.0, 4.5))
    ax2.plot(t, S1, label="S1 polar", lw=1.2)
    ax2.plot(t, S2, label="S2 director", lw=1.2)
    ax2.set_xlabel("t (s)")
    ax2.set_ylabel("order parameter")
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    fig2.tight_layout()
    fig2.savefig(png_path.with_name(png_path.stem + "_order.png"), dpi=160)
    plt.close(fig2)


# =============================================================================
# CLI
# =============================================================================



def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Research-grade inertial simulation of dipolar compass-needle arrays.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    grid_cfg = CFG.numerics.compass_engine.grid
    time_cfg = CFG.numerics.compass_engine.time
    tol_cfg = CFG.numerics.compass_engine.tolerances
    phys_cfg = CFG.physics.compass_engine
    run_cfg = CFG.run.compass_engine

    g = p.add_argument_group("lattice geometry")
    g.add_argument("--geometry", choices=["square", "triangular", "honeycomb"], default=grid_cfg.geometry)
    g.add_argument("--N", type=int, default=grid_cfg.N, help="number of rows or nominal honeycomb height")
    g.add_argument("--M", type=int, default=grid_cfg.M, help="number of columns or nominal honeycomb width")
    g.add_argument("--R", type=float, default=R_DEFAULT, help="half the centre-to-centre distance between pivots [m]; default gives 13 mm spacing")
    g.add_argument("--needle_frac", type=float, default=phys_cfg.needle_frac_legacy, help="legacy blade length fraction of 2R, used only with --use_legacy_size_from_R")
    g.add_argument("--needle_len", type=float, default=NEEDLE_LEN_DEFAULT, help="physical needle blade length [m]")
    g.add_argument("--needle_width", type=float, default=NEEDLE_WIDTH_DEFAULT, help="physical needle blade width [m]")
    g.add_argument("--use_legacy_size_from_R", type=int, choices=[0, 1], default=phys_cfg.use_legacy_size_from_R, help="if 1, set blade size from needle_frac*2R; if 0, use explicit needle_len/needle_width")

    g = p.add_argument_group("needle physical properties")
    g.add_argument("--moment", type=float, default=None, help="override magnetic moment [A m^2]")
    g.add_argument("--inertia", type=float, default=None, help="override moment of inertia [kg m^2]")
    g.add_argument("--needle_thickness", type=float, default=NEEDLE_THICKNESS_DEFAULT)
    g.add_argument("--steel_density", type=float, default=STEEL_DENSITY_DEFAULT)
    g.add_argument("--steel_Ms", type=float, default=STEEL_MS_SATURATION_DEFAULT, help="saturation magnetization [A/m]")
    g.add_argument("--steel_Bsat", type=float, default=None, help="optional Bsat [T], overrides steel_Ms via Ms=Bsat/mu0")
    g.add_argument("--pivot_radius", type=float, default=phys_cfg.pivot_radius)
    g.add_argument("--pivot_thickness", type=float, default=phys_cfg.pivot_thickness)
    g.add_argument("--pivot_density", type=float, default=phys_cfg.pivot_density)
    g.add_argument("--pivot_mass", type=float, default=None)
    g.add_argument("--damping", type=float, default=DAMPING_DEFAULT)
    g.add_argument("--damping_noise", type=float, default=phys_cfg.damping_noise, help="relative uniform random damping variation per needle")

    g = p.add_argument_group("time integration")
    g.add_argument("--t_sim", type=float, default=phys_cfg.t_sim, help="simulation time [s], except FORC and demag modes")
    g.add_argument("--dt_factor", type=float, default=time_cfg.dt_factor, help="dt/T0")
    g.add_argument("--noise", type=float, default=phys_cfg.noise, help="initial angular noise amplitude [rad]")
    g.add_argument("--seed", type=int, default=run_cfg.seed)
    g.add_argument("--log_every", type=int, default=time_cfg.log_every, help="write one CSV row every this many integration steps")
    g.add_argument("--flip_angle_deg", type=float, default=tol_cfg.flip_angle_deg, help="rest-angle displacement threshold for the committed flip_angle channel")
    g.add_argument("--flip_band_deg", type=float, default=tol_cfg.flip_band_deg, help="Schmitt dead-band half-width around the perpendicular to the drive axis [deg]")
    g.add_argument("--flip_dwell_T0", type=float, default=tol_cfg.flip_dwell_T0, help="dwell time required to commit a flip, in units of T0")
    g.add_argument("--flip_settle_frac", type=float, default=tol_cfg.flip_settle_frac, help="|omega| settling threshold as a fraction of omega0 for the flip_angle channel")
    g.add_argument("--event_log", action="store_true", default=run_cfg.event_log, help="write per-event CSV (step,t,needle_id,channel,theta) for offline avalanche clustering")
    g.add_argument("--dt_guard_alpha", type=float, default=tol_cfg.dt_guard_alpha, help="stability monitor threshold on max_i sqrt(m|B_i|/I)*dt")
    g.add_argument("--dt_guard_substep", action="store_true", default=run_cfg.dt_guard_substep, help="re-integrate flagged steps with 4 global sub-steps (breaks strict symplecticity; exploratory)")

    g = p.add_argument_group("dipolar cutoff and boundaries")
    g.add_argument("--cutoff_shells", type=float, default=grid_cfg.cutoff_shells, help="cutoff in units of r_nn")
    g.add_argument("--cutoff_m", type=float, default=None, help="absolute cutoff [m], overrides cutoff_shells")
    g.add_argument("--pbc", action="store_true", default=grid_cfg.pbc, help="periodic boundary conditions using finite images")
    g.add_argument("--n_images", type=int, default=grid_cfg.n_images, help="number of periodic images in each direction")
    g.add_argument("--tensor_mem_limit_gb", type=float, default=grid_cfg.tensor_mem_limit_gb)
    g.add_argument("--float32", action="store_true", default=grid_cfg.float32, help="use float32 tensor/state to save memory; float64 is preferred for final data")

    g = p.add_argument_group("field protocol")
    g.add_argument("--field_mode", choices=["static", "hysteresis", "forc", "sine", "step_up", "step_pos", "step_down", "step_neg", "pulse", "demag_rot", "demag_linear"], default=phys_cfg.field_mode)
    g.add_argument("--B_ext", type=float, default=None, help="field amplitude [T]; if omitted, B_max_factor*B_ref")
    g.add_argument("--B_max_factor", type=float, default=phys_cfg.B_max_factor, help="B_ext = factor * B_ref if B_ext omitted")
    g.add_argument("--phi_ext_deg", type=float, default=phys_cfg.phi_ext_deg, help="field direction")
    g.add_argument("--field_freq", type=float, default=phys_cfg.field_freq, help="sine frequency [Hz]")
    g.add_argument("--field_delay", type=float, default=phys_cfg.field_delay, help="delay for step/pulse protocols [s]")
    g.add_argument("--t_pulse", type=float, default=None, help="pulse duration [s]")
    g.add_argument("--hyst_spacing", choices=["linear", "log"], default=phys_cfg.hyst_spacing)
    g.add_argument("--hyst_log_k", type=float, default=phys_cfg.hyst_log_k)

    g = p.add_argument_group("FORC protocol")
    g.add_argument("--forc_Br_min", type=float, default=None, help="minimum reversal field [T]; default -Bmax")
    g.add_argument("--forc_n_curves", type=int, default=phys_cfg.forc.n_curves)
    g.add_argument("--forc_t_sat", type=float, default=phys_cfg.forc.t_sat)
    g.add_argument("--forc_t_ramp_down", type=float, default=phys_cfg.forc.t_ramp_down)
    g.add_argument("--forc_t_ramp_up", type=float, default=phys_cfg.forc.t_ramp_up)
    g.add_argument("--forc_rate", type=float, default=None, help="field ramp rate [T/s]; overrides fixed FORC ramp times")

    g = p.add_argument_group("demagnetization")
    g.add_argument("--demag_freq", type=float, default=phys_cfg.demag.freq)
    g.add_argument("--demag_cycles", type=int, default=phys_cfg.demag.cycles)
    g.add_argument("--t_relax_after", type=float, default=phys_cfg.demag.t_relax_after, help="extra relaxation after demag mode")

    g = p.add_argument_group("output and performance")
    g.add_argument("--out_dir", default=run_cfg.out_dir)
    g.add_argument("--tag", default=run_cfg.tag)
    g.add_argument("--use_gpu", action="store_true", default=run_cfg.use_gpu)
    g.add_argument("--progress", action="store_true", default=run_cfg.progress)
    g.add_argument("--verbose", action="store_true", default=run_cfg.verbose)
    g.add_argument("--make_plot", action="store_true", default=run_cfg.make_plot)
    g.add_argument("--png_dpi", type=int, default=run_cfg.png_dpi, help="DPI for automatic lattice PNG images")
    g.add_argument("--png_transparent", action="store_true", default=run_cfg.png_transparent, help="save lattice PNGs with transparent background")
    g.add_argument("--png_with_axes", action="store_true", default=run_cfg.png_with_axes, help="include axes/grid/title in lattice PNGs")
    g.add_argument("--png_no_panel_titles", action="store_true", default=run_cfg.png_no_panel_titles, help="remove Initial/Final panel titles in side-by-side PNG")
    g.add_argument("--domain_tol_deg", type=float, default=tol_cfg.domain_tol_deg)

    return p



def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    csv_path, meta_path, state_path = run_simulation(args)
    print(f"CSV      : {csv_path}")
    print(f"Metadata : {meta_path}")
    print(f"Final NPZ: {state_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
