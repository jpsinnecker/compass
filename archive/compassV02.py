#!/usr/bin/env python3
"""
compassV02.py — research-grade simulation of dipolar compass-needle arrays
Version: V2.0.2, July 2026

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

3. Avalanche activity is exported directly:
   flip_field counts sign changes relative to the selected drive axis.
   flip_angle counts rotations larger than --flip_angle_deg in one time step.
   These are accumulated over each logging interval and should be used for
   avalanche statistics instead of the continuous order parameter S1.

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


# =============================================================================
# Physical constants and defaults
# =============================================================================

MU0_OVER_4PI = 1.0e-7  # T m / A
STEEL_DENSITY_DEFAULT = 7850.0  # kg/m^3
STEEL_MS_SATURATION_DEFAULT = 1.59e6  # A/m, approximately Bsat=2.0 T / mu0
# Default physical apparatus geometry used in the project report.
# The lattice parameter is defined by the centre-to-centre spacing between
# neighbouring pivots. In the code the historical variable R is kept as half
# this distance because the lattice generators use d = 2R.
CENTER_DISTANCE_DEFAULT = 0.013  # m, centre-to-centre pivot distance
R_DEFAULT = 0.5 * CENTER_DISTANCE_DEFAULT  # m
NEEDLE_LEN_DEFAULT = 0.010  # m, physical blade length
NEEDLE_WIDTH_DEFAULT = 0.003  # m, physical blade width
NEEDLE_THICKNESS_DEFAULT = 0.0004  # m
DAMPING_DEFAULT = 5.0e-8  # N m s / rad


# =============================================================================
# Small numerical helpers
# =============================================================================


def wrap_angle(a):
    """Wrap angle(s) to (-pi, pi]. Works for NumPy or CuPy arrays."""
    return (a + math.pi) % (2.0 * math.pi) - math.pi


def wrap_angle_np(a: np.ndarray) -> np.ndarray:
    return (a + np.pi) % (2.0 * np.pi) - np.pi


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
        needle_width_eff = 0.22 * needle_len_eff
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
        dtheta = np.abs(wrap_angle_np(theta_cpu[i] - theta_cpu[neigh]))
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
    domain_tol_deg: float



def run_simulation(args) -> Tuple[Path, Path, Path]:
    out_dir = ensure_dir(args.out_dir)
    data_dir = ensure_dir(out_dir / "data")
    meta_dir = ensure_dir(out_dir / "meta")
    state_dir = ensure_dir(out_dir / "states")

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

    # Initial field and torque.
    f0 = protocol.at(0.0)
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    def torque_and_field(th, bx, by):
        mx = moment * xp.cos(th)
        my = moment * xp.sin(th)
        Bdx = tensor.Axx @ mx + tensor.Axy @ my
        Bdy = tensor.Axy @ mx + tensor.Ayy @ my
        return mx * (Bdy + by) - my * (Bdx + bx)

    tau = torque_and_field(theta, f0.bx, f0.by)
    sigma_prev = xp.where(xp.cos(theta - phi) >= 0.0, 1, -1)
    flip_field_acc = 0
    flip_angle_acc = 0

    invI = 1.0 / inertia
    dt2_half = 0.5 * dt * dt
    b_half = damping * dt * 0.5 * invI
    coeff_omega = (1.0 - b_half) / (1.0 + b_half)
    coeff_tau = dt * 0.5 * invI / (1.0 + b_half)

    tag = args.tag or f"{args.geometry}_{args.field_mode}_N{args.N}_M{args.M}_seed{seed}"
    csv_path = data_dir / f"{tag}.csv"
    meta_path = meta_dir / f"{tag}.json"
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
        domain_tol_deg=args.domain_tol_deg,
    )

    metadata = {
        "program": "compassV02.py",
        "version": "2.0.2",
        "created_unix_time": time.time(),
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
        },
        "notes": {
            "S1": "polar order |<exp(i theta)>|",
            "S2": "nematic/director order |<exp(2 i theta)>|",
            "flip_field": "accumulated sign changes relative to the drive axis since previous log row",
            "flip_angle": "accumulated rotations larger than flip_angle_deg since previous log row",
            "E_dip": "-0.5 sum_i m_i dot B_i^dip",
            "E_ext": "-sum_i m_i dot B_ext",
            "E_kin": "0.5 I sum_i omega_i^2",
        },
    }

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

        for step in range(1, n_steps + 1):
            t = step * dt
            f = protocol.at(t)

            theta_old = theta
            accel = (tau - damping * omega) * invI
            theta_new = theta + omega * dt + accel * dt2_half
            theta_new = (theta_new + math.pi) % (2.0 * math.pi) - math.pi

            tau_new = torque_and_field(theta_new, f.bx, f.by)
            omega_new = omega * coeff_omega + (tau + tau_new) * coeff_tau

            sigma_new = xp.where(xp.cos(theta_new - phi) >= 0.0, 1, -1)
            flip_field_acc += int(safe_float(xp.sum(sigma_new != sigma_prev)))
            dtheta_step = (theta_new - theta_old + math.pi) % (2.0 * math.pi) - math.pi
            flip_angle_acc += int(safe_float(xp.sum(xp.abs(dtheta_step) >= flip_threshold)))
            sigma_prev = sigma_new

            theta = theta_new
            omega = omega_new
            tau = tau_new

            if step % log_every == 0 or step == n_steps:
                write_row(step, t, f)

            if args.progress and step % max(1, n_steps // 100) == 0:
                elapsed = time.perf_counter() - t_start
                rate = step / max(elapsed, 1e-12)
                print(f"\r{100*step/n_steps:6.2f}%  step {step}/{n_steps}  {rate:8.0f} steps/s", end="", flush=True)

    if args.progress:
        print()
    backend.sync()

    theta_cpu = backend.to_cpu(theta)
    omega_cpu = backend.to_cpu(omega)
    domain_stats = domain_statistics(theta_cpu, geom.xs, geom.ys, geom.r_nn, args.domain_tol_deg)
    metadata["domain_statistics_final"] = domain_stats
    metadata["runtime_s"] = time.perf_counter() - t_start

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

    if args.make_plot:
        try:
            make_quick_plot(csv_path, out_dir / f"{tag}_quicklook.png")
        except Exception as exc:
            print(f"[warning] quick plot failed: {exc}")

    return csv_path, meta_path, state_path


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

    g = p.add_argument_group("lattice geometry")
    g.add_argument("--geometry", choices=["square", "triangular", "honeycomb"], default="square")
    g.add_argument("--N", type=int, default=16, help="number of rows or nominal honeycomb height")
    g.add_argument("--M", type=int, default=16, help="number of columns or nominal honeycomb width")
    g.add_argument("--R", type=float, default=R_DEFAULT, help="half the centre-to-centre distance between pivots [m]; default gives 13 mm spacing")
    g.add_argument("--needle_frac", type=float, default=0.80, help="legacy blade length fraction of 2R, used only with --use_legacy_size_from_R")
    g.add_argument("--needle_len", type=float, default=NEEDLE_LEN_DEFAULT, help="physical needle blade length [m]")
    g.add_argument("--needle_width", type=float, default=NEEDLE_WIDTH_DEFAULT, help="physical needle blade width [m]")
    g.add_argument("--use_legacy_size_from_R", type=int, choices=[0, 1], default=0, help="if 1, set blade size from needle_frac*2R; if 0, use explicit needle_len/needle_width")

    g = p.add_argument_group("needle physical properties")
    g.add_argument("--moment", type=float, default=None, help="override magnetic moment [A m^2]")
    g.add_argument("--inertia", type=float, default=None, help="override moment of inertia [kg m^2]")
    g.add_argument("--needle_thickness", type=float, default=NEEDLE_THICKNESS_DEFAULT)
    g.add_argument("--steel_density", type=float, default=STEEL_DENSITY_DEFAULT)
    g.add_argument("--steel_Ms", type=float, default=STEEL_MS_SATURATION_DEFAULT, help="saturation magnetization [A/m]")
    g.add_argument("--steel_Bsat", type=float, default=None, help="optional Bsat [T], overrides steel_Ms via Ms=Bsat/mu0")
    g.add_argument("--pivot_radius", type=float, default=0.0)
    g.add_argument("--pivot_thickness", type=float, default=0.0)
    g.add_argument("--pivot_density", type=float, default=8500.0)
    g.add_argument("--pivot_mass", type=float, default=None)
    g.add_argument("--damping", type=float, default=DAMPING_DEFAULT)
    g.add_argument("--damping_noise", type=float, default=0.0, help="relative uniform random damping variation per needle")

    g = p.add_argument_group("time integration")
    g.add_argument("--t_sim", type=float, default=2.0, help="simulation time [s], except FORC and demag modes")
    g.add_argument("--dt_factor", type=float, default=0.04, help="dt/T0")
    g.add_argument("--noise", type=float, default=1.5, help="initial angular noise amplitude [rad]")
    g.add_argument("--seed", type=int, default=None)
    g.add_argument("--log_every", type=int, default=10, help="write one CSV row every this many integration steps")
    g.add_argument("--flip_angle_deg", type=float, default=90.0, help="threshold for large-rotation activity counter")

    g = p.add_argument_group("dipolar cutoff and boundaries")
    g.add_argument("--cutoff_shells", type=float, default=3.5, help="cutoff in units of r_nn")
    g.add_argument("--cutoff_m", type=float, default=None, help="absolute cutoff [m], overrides cutoff_shells")
    g.add_argument("--pbc", action="store_true", help="periodic boundary conditions using finite images")
    g.add_argument("--n_images", type=int, default=1, help="number of periodic images in each direction")
    g.add_argument("--tensor_mem_limit_gb", type=float, default=6.0)
    g.add_argument("--float32", action="store_true", help="use float32 tensor/state to save memory; float64 is preferred for final data")

    g = p.add_argument_group("field protocol")
    g.add_argument("--field_mode", choices=["static", "hysteresis", "forc", "sine", "step_up", "step_pos", "step_down", "step_neg", "pulse", "demag_rot", "demag_linear"], default="static")
    g.add_argument("--B_ext", type=float, default=None, help="field amplitude [T]; if omitted, B_max_factor*B_ref")
    g.add_argument("--B_max_factor", type=float, default=8.0, help="B_ext = factor * B_ref if B_ext omitted")
    g.add_argument("--phi_ext_deg", type=float, default=0.0, help="field direction")
    g.add_argument("--field_freq", type=float, default=1.0, help="sine frequency [Hz]")
    g.add_argument("--field_delay", type=float, default=0.0, help="delay for step/pulse protocols [s]")
    g.add_argument("--t_pulse", type=float, default=None, help="pulse duration [s]")
    g.add_argument("--hyst_spacing", choices=["linear", "log"], default="linear")
    g.add_argument("--hyst_log_k", type=float, default=5.0)

    g = p.add_argument_group("FORC protocol")
    g.add_argument("--forc_Br_min", type=float, default=None, help="minimum reversal field [T]; default -Bmax")
    g.add_argument("--forc_n_curves", type=int, default=30)
    g.add_argument("--forc_t_sat", type=float, default=0.05)
    g.add_argument("--forc_t_ramp_down", type=float, default=0.10)
    g.add_argument("--forc_t_ramp_up", type=float, default=0.20)
    g.add_argument("--forc_rate", type=float, default=None, help="field ramp rate [T/s]; overrides fixed FORC ramp times")

    g = p.add_argument_group("demagnetization")
    g.add_argument("--demag_freq", type=float, default=2.0)
    g.add_argument("--demag_cycles", type=int, default=20)
    g.add_argument("--t_relax_after", type=float, default=2.0, help="extra relaxation after demag mode")

    g = p.add_argument_group("output and performance")
    g.add_argument("--out_dir", default="compassV2_output")
    g.add_argument("--tag", default=None)
    g.add_argument("--use_gpu", action="store_true")
    g.add_argument("--progress", action="store_true")
    g.add_argument("--verbose", action="store_true")
    g.add_argument("--make_plot", action="store_true")
    g.add_argument("--domain_tol_deg", type=float, default=15.0)

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
