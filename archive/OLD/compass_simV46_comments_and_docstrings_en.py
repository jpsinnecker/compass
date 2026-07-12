"""
============================================================
compass_sim.py — Simulation of a compass-needle lattice
============================================================

Models a 2D grid of classical magnetic dipoles, represented as compass
needles, that interact through the magnetic field each needle produces at
its neighbors. The dynamics are inertial, using Newton's second law for
rotation with no pivot friction and viscous air damping. The integrator is
Velocity-Verlet. All physical quantities use SI units.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
COMMAND-LINE PARAMETERS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

┌─ LATTICE GEOMETRY ────────────────────────────────────────┐
│                                                            │
│  --geometry  square | triangular | honeycomb               │
│              Lattice type. Default: square                 │
│              · square     : rectangular grid               │
│              · triangular : equilateral triangular grid    │
│              · honeycomb  : honeycomb, with hexagonal holes│
│                                                            │
│  --N         int   Number of needle ROWS. Default: 8       │
│  --M         int   Number of needle COLUMNS. Default: 8    │
│                                                            │
│  --R         float Radius of each needle circle [m].       │
│                    Neighbor distance = 2R.                 │
│                    Default: 0.025 m (2.5 cm)               │
│                                                            │
│  --needle_frac float Needle length as a fraction of        │
│                    the diameter 2R (0.0-1.0).              │
│                    Default: 0.80 -> needle = 80% of 2R     │
└────────────────────────────────────────────────────────────┘

┌─ PHYSICS AND SIMULATION ──────────────────────────────────┐
│                                                            │
│  --moment    float Magnetic moment of each needle [A·m²]  │
│                    Default: 0.1 (table compass, about 5 cm)│
│                    Ref: pocket ≈ 0.01 | nautical ≈ 1.0     │
│                                                            │
│  --inertia   float Moment of inertia [kg·m²]. If omitted, │
│                    it is computed automatically from the   │
│                    geometry: steel sheet using R,          │
│                    needle_frac, --needle_thickness,        │
│                    and --steel_density.                    │
│                                                            │
│  --damping   float Viscous air damping [N·m·s/rad]        │
│                    Controls the quality factor Q:          │
│                    Q = omega_0·I/b (high Q = more          │
│                    oscillatory motion). Default: 5e-8      │
│                    (Q≈25, realistic compass). For a smooth │
│                    B_ext=0.1 T run, use 8e-6 (Q≈4).        │
│                                                            │
│  --t_sim     float Total physical simulation time [s].     │
│                    Sum of all integrated dt steps,          │
│                    equivalent to what a real stopwatch     │
│                    would measure while observing the       │
│                    moving needles. The video displays this │
│                    physical time. The simulation stops     │
│                    earlier if S=1.00 or at rest.           │
│                    Default: 2.0 s                          │
│                                                            │
│  --dt_factor float Fraction of the natural period T0 used  │
│                    as the integration step (0.02-0.10).    │
│                    Smaller = more accurate, slower.        │
│                    Default: 0.05                           │
│                                                            │
│  --noise     float Initial-noise amplitude [rad].          │
│                    0 = all needles point toward +x         │
│                    pi ≈ 3.14 = fully random orientation    │
│                    Default: 1.5                            │
│                                                            │
│  --seed      int   Random-number-generator seed.           │
│                    Ensures reproducibility.                │
│                    Default: 42                             │
└────────────────────────────────────────────────────────────┘

┌─ UNIFORM EXTERNAL FIELD (SI units) ────────────────────────┐
│  Two ways to specify the field. Do not use both together.  │
│                                                            │
│  Form A — magnitude + angle, recommended:                  │
│  --B_ext     float Field magnitude [T].                    │
│                    0.0      = no field (default)           │
│                    50e-6    = Earth's field (≈50 µT)       │
│                    1e-3     = fridge magnet at 5 cm        │
│                    0.1      = strong field (aligns all)    │
│                                                            │
│  --phi_ext   float Field direction [degrees].              │
│                    0   = right (+x), default               │
│                    90  = up (+y)                           │
│                    180 = left (-x)                         │
│                    270 = down (-y)                         │
│                    Counterclockwise sense. Default: 0.0    │
│                                                            │
│  Form B — Cartesian components, overrides A:               │
│  --ext_Bx    float x component of the field [T]            │
│  --ext_By    float y component of the field [T]            │
└────────────────────────────────────────────────────────────┘

┌─ OUTPUT ──────────────────────────────────────────────────┐
│  PNG files are always generated in the current directory:  │
│    compass_initial.png      initial state                  │
│    compass_equilibrium.png  final state                    │
│    compass_comparison.png   side-by-side comparison        │
│    compass_order_param.png  order parameter S(t)           │
│                                                            │
│  --video     str   Path of the MP4 video to generate.      │
│                    Requires ffmpeg. If the file already    │
│                    exists, saves as name0001.mp4,          │
│                    name0002.mp4, etc.                     │
│                    Example: --video simulation.mp4         │
│                                                            │
│  --frame_every int Save one frame every N steps.           │
│                    Smaller = smoother video, slower run.   │
│                    Default: 5                              │
│                                                            │
│  --fps       int   Frames per second of the MP4 video.     │
│                    Default: 24                             │
│                                                            │
│  --keep_frames     If present, keeps the folder of         │
│                    intermediate PNG frames after MP4       │
│                    generation.                             │
└────────────────────────────────────────────────────────────┘

┌─ CONTROLS DURING THE SIMULATION ──────────────────────────┐
│  The simulation stops automatically when:                  │
│    · S = 1.00, all needles aligned                         │
│    · The lattice is at rest, omega_max -> 0                │
│    · The requested t_sim time is reached                   │
│    · Ctrl+I (Tab) is pressed in the terminal               │
│      -> interrupts and saves the video immediately         │
└────────────────────────────────────────────────────────────┘

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
EXAMPLES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  # Defaults — 8x8 square lattice with no field
  python compass_sim.py

  # Honeycomb with Earth's field
  python compass_sim.py --geometry honeycomb --N 10 --M 10 --B_ext 50e-6

  # Triangular lattice with a 0.1 T field at 45°, smooth motion
  python compass_sim.py --geometry triangular --N 10 --M 10 \
      --B_ext 0.1 --phi_ext 45 --damping 8e-6 --t_sim 2.0

  # Video with larger needles and many oscillations
  python compass_sim.py --R 0.03 --needle_frac 0.85 --damping 1e-9 \
      --t_sim 5.0 --frame_every 2 --fps 30 --video sim.mp4

  # Field through Cartesian components
  python compass_sim.py --ext_Bx 0.05 --ext_By -0.05

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Dependencies: numpy, matplotlib, ffmpeg for video
Optional GPU acceleration: cupy-cuda12x, or a version compatible with your CUDA
  pip install cupy-cuda12x[ctk]
  The [ctk] suffix also installs CUDA Toolkit headers, which are required
  to compile kernels at runtime. Without them, the program falls back to CPU
  with the error "Failed to find CUDA headers".
  The backend is detected automatically. If CuPy is absent or incomplete,
  the program uses CPU/NumPy.
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
from matplotlib.collections import PatchCollection, LineCollection
import argparse
import time
import sys

# ── GPU acceleration backend (optional) ──────────────────────────────────
# If the package 'cupy' is installed AND an NVIDIA GPU (CUDA) estiver
# available, uses CuPy as a NumPy replacement for heavy computations
# (dipolar field). CuPy has an API that is almost identical to NumPy, so the same
# vectorized code runs on the GPU without logic changes.
# If cupy is not installed, there is no GPU, or the CUDA headers
# Toolkit are not present (needed to compile kernels in
# runtime), automatically falls back to NumPy (CPU) — the programa
# works on any machine, even with cupy "partially" installed.
_GPU_AVAILABLE = False
_GPU_ERROR_MSG = None
_xp = np   # _xp = "array module" active: np (CPU) or cp (GPU)
try:
    import cupy as cp
    # Test 1: is there the visible CUDA device?
    cp.cuda.Device(0).compute_capability
    # Test 2 (critical): actually tries to COMPILE and execute the simple kernel.
    # cp.cuda.Device(0).compute_capability only queries the driver — does not detect
    # missing headers of the CUDA Toolkit, which are required for the NVRTC
    # compilar kernels elementwise (cos, sin, etc.) in runtime.
    # Without this real test, the error only appears in the middle of the simulation.
    _test_arr = cp.array([0.0, 1.0])
    _ = cp.cos(_test_arr)   # forces JIT compilation of the real kernel
    cp.cuda.Stream.null.synchronize()
    _xp = cp
    _GPU_AVAILABLE = True
except Exception as _e:
    _xp = np
    _GPU_AVAILABLE = False
    _GPU_ERROR_MSG = str(_e).strip().splitlines()[-1] if str(_e).strip() else type(_e).__name__


def _to_cpu(arr):
    """Converts an array from the active backend, GPU or CPU, to a plain NumPy array."""
    if _GPU_AVAILABLE and hasattr(arr, 'get'):
        return arr.get()
    return arr


def _to_backend(arr):
    """Converts a NumPy array to the active backend, GPU or CPU."""
    if _GPU_AVAILABLE:
        return cp.asarray(arr)
    return np.asarray(arr)


# ── configure output for immediate flushing and correct line endings ──────────────
# On macOS with zsh, keyboard raw mode (used by Ctrl+I) may disable
# the processing of \n → \r\n in the terminal. Forcing line_buffering and flush
# ensures that each print() appears correctly on its own line.
sys.stdout.reconfigure(line_buffering=True)

def _print(*args, **kwargs):
    """Replacement for print() that guarantees immediate flushing and correct CRLF line endings, even when the terminal is in raw mode, as on macOS/zsh with Ctrl+I enabled."""
    import sys
    msg = ' '.join(str(a) for a in args)
    sys.stdout.write(msg + '\r\n')
    sys.stdout.flush()


def _progress_ansi_ok():
    """Returns True if stdout is an interactive terminal (TTY). This enables the live two-line panel using ANSI codes. For redirected output, such as logs or pipes, it returns False and the code uses the legacy behavior without escape codes that would pollute the file."""
    import sys
    if _PROGRESS_REGION['ansi'] is None:
        try:
            _PROGRESS_REGION['ansi'] = bool(sys.stdout.isatty())
        except Exception:
            _PROGRESS_REGION['ansi'] = False
    return _PROGRESS_REGION['ansi']


# state of the live progress-bar panel (V44):
#   'open' : True while the two lines (bar + status) are on the screen
#   'ansi' : TTY detection cache (None = not yet checked)
_PROGRESS_REGION = {'open': False, 'ansi': None}


def _print_progress_bar(frac, prefix="", suffix="", status_line=None,
                        bar_width=40):
    """Prints or updates the progress panel (V44): two live lines rewritten in place at each update.

Example:
    Integrating [CPU] [########            ]  42.3%  step 970/2294
    [step 970/2294]  t=2.4310s  B=+0.0000mT  S=0.4415  w_max=3.1rad/s

In an interactive terminal (TTY), ANSI codes move the cursor and clear lines so the two lines are redrawn in place instead of accumulating new lines. For redirected output (non-TTY), the function uses the legacy behavior: the bar is updated with '\r' on a single line, and the caller prints status as normal log lines.

Parameters
----------
frac        : completed fraction, between 0.0 and 1.0
prefix      : text before the bar, e.g. "Integrating [CPU] "
suffix      : text after the percentage, e.g. "step 524/1007"
status_line : second panel line, or None
bar_width   : number of characters in the bar itself

Note: any normal _print() message emitted while the panel is open should come after _print_progress_bar_finish(), which moves the cursor below the panel and marks it as closed."""
    import sys
    frac = max(0.0, min(1.0, frac))
    n_filled = int(round(frac * bar_width))
    bar = '#' * n_filled + ' ' * (bar_width - n_filled)
    pct = frac * 100.0
    bar_line = f"{prefix}[{bar}] {pct:5.1f}%  {suffix}"

    if _progress_ansi_ok():
        status = status_line if status_line is not None else ""
        # Protocol with cursor invariant (V44): at the END of each drawing,
        # the cursor remains parked at COLUMN 0 OF THE BAR LINE (top of the
        # panel). Each update then redraws downward:
        #   clear line (\x1b[K) + bar + moves down ("\n\r", robust with or
        #   without translation ONLCR by the terminal driver) + clears + status +
        #   moves up (\x1b[1THE) + returns to column 0 (\r).
        # This is independent of the position where the text ends (without depending
        # on \n returning to column 0, which varies between terminals/ptys) and
        # works at the bottom of the screen (the "\n" causes the required scroll).
        seq = ("\x1b[K" + bar_line + "\n\r" +
               "\x1b[K" + status + "\x1b[1A\r")
        sys.stdout.write(seq)
        sys.stdout.flush()
        _PROGRESS_REGION['open'] = True
    else:
        # legacy (not-TTY): bar overwriting its own line via \r
        sys.stdout.write("\r" + bar_line)
        sys.stdout.flush()


def _print_progress_bar_finish():
    """Closes the progress panel: moves the cursor to a new clean line below the panel, keeps the final two lines visible, and marks the panel as closed. In ANSI mode it is safe to call this when the panel is already closed, because it emits nothing and avoids spurious blank lines."""
    import sys
    if _progress_ansi_ok():
        if _PROGRESS_REGION['open']:
            # cursor is at column 0 of the bar line: move down two lines
            # (bar -> status -> new line below the panel)
            sys.stdout.write('\n\r\n\r')
            sys.stdout.flush()
            _PROGRESS_REGION['open'] = False
    else:
        sys.stdout.write('\r\n')
        sys.stdout.flush()

# ══════════════════════════════════════════════════════════════════════════════
# 1. PHYSICAL CONSTANTS (International System — SI)
# ══════════════════════════════════════════════════════════════════════════════

# μ₀/(4π) = 1×10⁻⁷  T·m/A  (vacuum permeability / 4π)
MU0_OVER_4PI = 1.0e-7   # T·m/A

# ── Default physical parameters of a table compass needle (~5 cm) ──────
#
# Magnetic moment
MOMENT_DEFAULT = 0.1      # A·m²   (typical table compass needle)
#
# Needle mass and geometry (thin steel bar, length L, circular cross section)
#   mass   m ≈ 0.5 g  = 5×10⁻⁴ kg
#   length L ≈ 5 cm = 0.05 m
#   Moment of inertia of thin bar about its center:
#       I = (1/12) · m · L²  =  (1/12) · 5×10⁻⁴ · (0.05)²  ≈ 1.04×10⁻⁷  kg·m²
INERTIA_DEFAULT = 1.0e-7  # kg·m²

# Density of common magnetic steel (carbon steel / silicon steel for use
# eletromagnetic, the used in transformer laminations or compasss
# of mesa). Typically ranges between 7700-7900 kg/m³; we use the value
# representative of the ordinary carbon steel.
STEEL_DENSITY_DEFAULT = 7850.0   # kg/m³
# Default thickness of a thin sheet cut from ordinary steel sheet
NEEDLE_THICKNESS_DEFAULT = 0.4e-3   # m  (0.4 mm)

# Saturation magnetization of a common magnetic steel (carbon steel / ferro
# doce). Corresponds to a saturatestion flux density Bsat ≈ 2.0 T,
# well-established value in the literature for iron alloys of alta
# permeability (typical range: 1.6-2.2 T; we use the value central
# representative). Relation: Bsat = μ0 · Ms  →  Ms = Bsat / μ0.
STEEL_MS_SATURATION_DEFAULT = 1.59e6   # A/m  (≈ Bsat of 2.0 T)


def compute_inertia_from_geometry(needle_len, needle_width, thickness,
                                  density=STEEL_DENSITY_DEFAULT):
    """Computes the moment of inertia of a needle modeled as a thin, homogeneous rectangular sheet rotating about the axis perpendicular to the sheet plane and passing through the center of mass. This is the z axis and represents the real rotation axis of a compass needle supported by a pivot.

Physical model
--------------
Rectangular sheet with length L, width w, thickness t, mass m, and density rho:

    m = rho · L · w · t                         [kg]

Moment of inertia of a thin rectangular plate about the axis perpendicular to its plane through the center:

    I = (1/12) · m · (L² + w²)                  [kg·m²]

This formula generalizes the thin-bar case, I = mL²/12, valid when w -> 0. Here the w² term is kept explicitly because real needles cut from sheet metal have non-negligible width. The thickness t enters only through the mass. A thin sheet has a negligible intrinsic inertia term in the direction perpendicular to its plane, m·t²/12, when t << L,w.

Parameters
----------
needle_len   : needle length L [m]
needle_width : needle width w [m]
thickness    : sheet thickness t [m]
density      : material density [kg/m³], default: common magnetic steel, 7850 kg/m³

Returns
-------
I : moment of inertia [kg·m²]"""
    mass = density * needle_len * needle_width * thickness   # [kg]
    I = (1.0 / 12.0) * mass * (needle_len**2 + needle_width**2)   # [kg·m²]
    return I


def compute_moment_from_geometry(needle_len, needle_width, thickness,
                                 Ms=STEEL_MS_SATURATION_DEFAULT):
    """Computes the magnetic moment of a needle modeled as a thin rectangular sheet of magnetic steel, saturated along its long axis, the north direction of the needle. In other words, the whole needle is assumed to be magnetized up to the material's physical limit, with all magnetic domains aligned along its length.

Physical model
--------------
For a saturated ferromagnetic material, the total magnetic moment is the product of saturation magnetization Ms and volume V:

    m = Ms · V                                  [A·m²]
    V = L · w · t                               [m³]

This is the same rectangular-sheet geometry L×w×t used in the inertia calculation. Here the full needle volume is assumed to contribute to magnetization, an ideal case of uniform saturation. Real needles may have magnetization slightly below the theoretical saturation because of shape demagnetization, but full saturation is the physical upper limit and a reasonable reference for a well-magnetized needle.

Parameters
----------
needle_len   : needle length L [m]
needle_width : needle width w [m]
thickness    : sheet thickness t [m]
Ms           : material saturation magnetization [A/m], default: common magnetic steel, about 1.59e6 A/m, corresponding to Bsat ≈ 2.0 T

Returns
-------
m : magnetic moment [A·m²]"""
    volume = needle_len * needle_width * thickness   # [m³]
    m = Ms * volume                                   # [A·m²]
    return m

#
# Viscous damping (air resistance to needle rotation).
# The effective Q depends of the dominant field B_eff = max(B_dipolar, B_externo):
#   Q = omega_0·I/b   where  omega_0 = sqrt(m·B_eff/I)
#
# With B_ext = 0.1 T:   omega_0 ≈ 316 rad/s  →  b = 8and-6 gives Q ≈ 4  (smooth)
# With B_ext = 0  :     omega_0 ≈  13 rad/s  →  b = 5and-8 gives Q ≈ 25 (compass)
#
# The default 5e-8 is suitable for simulating with no external field or weak field.
# For strong fields (> 1 mT), use --damping 1and-6 the 1and-5 for smooth motion.
DAMPING_DEFAULT = 5.0e-8  # N·m·s/rad


# ══════════════════════════════════════════════════════════════════════════════
# 2. MAGNETIC DIPOLAR FIELD (SI)
# ══════════════════════════════════════════════════════════════════════════════

def dipole_field_2d(rx, ry, theta_src, moment):
    """Computes the magnetic field produced by a dipole in the XY plane, using SI units.

Exact 3D magnetic-dipole field, evaluated in the z=0 plane:

    B = (mu0/4pi) · [ 3(m_hat·r_hat)r_hat - m_hat ] / r³

In Cartesian components, with m = moment·(cos theta, sin theta, 0):

    Bx = (mu0/4pi) · [ 3(m·r)·rx / r⁵ - mx / r³ ]
    By = (mu0/4pi) · [ 3(m·r)·ry / r⁵ - my / r³ ]

SI units:
    moment [A·m²]
    rx, ry [m]
    Bx, By [T]

Parameters
----------
rx, ry    : components of r = receiver point - source dipole [m]
theta_src : angle of the source dipole relative to the +x axis [rad]
moment    : magnitude of the magnetic dipole moment of the source needle [A·m²]

Returns
-------
(Bx, By) : magnetic-field components at (rx, ry) [T]"""
    r2 = rx**2 + ry**2
    if r2 < 1e-24:
        # avoids division by zero when the point coincides with the source
        return 0.0, 0.0

    r  = np.sqrt(r2)          # distance  [m]
    r5 = r2 * r2 * r          # r⁵  [m⁵]

    mx = moment * np.cos(theta_src)   # x component of the moment  [THE·m²]
    my = moment * np.sin(theta_src)   # y component of the moment  [THE·m²]

    mdotr = mx * rx + my * ry         # m · r  [THE·m³]

    Bx = MU0_OVER_4PI * (3.0 * mdotr * rx / r5  -  mx / (r2 * r))
    By = MU0_OVER_4PI * (3.0 * mdotr * ry / r5  -  my / (r2 * r))
    return Bx, By   # [T]


# ══════════════════════════════════════════════════════════════════════════════
# 3. TOTAL FIELD ON THE NEEDLE (SI)
# ══════════════════════════════════════════════════════════════════════════════

def total_field_on(i, j, thetas, xs, ys, cutoff, moment,
                   pbc=False, Lx=None, Ly=None, n_images=1):
    """Sums the dipolar field from all neighboring needles onto needle (i, j).

For computational efficiency, only needles within the cutoff radius are considered. This is physically justified because the dipolar field decays as 1/r³, so distant contributions are negligible.

Periodic boundary conditions (PBC)
----------------------------------
When pbc=True, the lattice is treated as a unit cell repeated infinitely in x and y, creating a periodic structure without edges. For each needle pair (i,j)-(ni,nj), the field is summed over all periodic replicas within ±n_images cells in x and y, not only the nearest image. With n_images=1, the default, each neighbor contributes through (2·1+1)² = 9 replicas: the original cell plus one replica on each side in x and y, equivalent to the classical minimum-image convention. Increasing n_images includes more distant replicas and is useful when the cutoff is large relative to the lattice period.

This removes edge effects: needles at the boundary see neighbors from the opposite side, and their more distant replicas, as if the structure were infinite.

Parameters
----------
i, j     : indices of the receiver needle in the N×M grid
thetas   : 2D array with current angles of all needles [rad]
xs, ys   : 2D arrays with fixed needle positions [m]
cutoff   : maximum interaction radius [m]
moment   : magnetic moment of each needle [A·m²]
pbc      : if True, applies periodic boundary conditions [bool]
Lx, Ly   : lattice periods in x and y [m], required only when pbc=True and computed in make_grid
n_images : number of periodic replicas to sum on each side in each direction, default 1. Total replicas considered per neighbor = (2·n_images+1)² cells.

Returns
-------
(Bx_tot, By_tot) : total field at needle (i, j) [T]"""
    Bx_tot, By_tot = 0.0, 0.0
    N, M = thetas.shape
    xi, yi = xs[i, j], ys[i, j]

    # periodic replica shifts to consider (0 = original cell)
    if pbc and Lx and Ly:
        img_range = range(-n_images, n_images + 1)
        x_shifts = [k * Lx for k in img_range]
        y_shifts = [k * Ly for k in img_range]
    else:
        x_shifts = [0.0]
        y_shifts = [0.0]

    for ni in range(N):
        for nj in range(M):
            if ni == i and nj == j and not pbc:
                continue  # without PBC: the needle does not interact with itself

            rx0 = xi - xs[ni, nj]   # [m]
            ry0 = yi - ys[ni, nj]   # [m]

            for dx_img in x_shifts:
                for dy_img in y_shifts:
                    # skips the needle itself in the original cell (zero distance)
                    if ni == i and nj == j and dx_img == 0.0 and dy_img == 0.0:
                        continue

                    rx = rx0 + dx_img
                    ry = ry0 + dy_img
                    dist = np.sqrt(rx*rx + ry*ry)

                    if dist > cutoff:
                        continue  # outside the cutoff radius: ignore

                    bx, by = dipole_field_2d(rx, ry, thetas[ni, nj], moment)
                    Bx_tot += bx
                    By_tot += by

    return Bx_tot, By_tot   # [T]

    return Bx_tot, By_tot   # [T]


# ══════════════════════════════════════════════════════════════════════════════
# 4. INERTIAL DYNAMICS (Newton's second law for rotation — no pivot friction)
# ══════════════════════════════════════════════════════════════════════════════

def _plot_hysteresis(log):
    """Plots the hysteresis curve M_proj(B) and saves it as PNG and CSV.

The plot shows the magnetization projected along the field direction as a function of the applied-field magnitude, the classical M-H hysteresis curve.

Parameters
----------
log : list of tuples (t, B_scalar, M_proj, S), generated during integration with field_mode='hysteresis' """
    log = np.array(log)
    t_arr  = log[:, 0]
    B_arr  = log[:, 1] * 1e3    # converts T → mT for readability
    M_arr  = log[:, 2]          # projected magnetization (dimensionless, ∈ [−1,1])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), facecolor='#1A1A2E')

    # ── left panel: M×B (hysteresis curve) ─────────────────────────
    ax1 = axes[0]
    ax1.set_facecolor('#0D1B2A')
    ax1.plot(B_arr, M_arr, color='#E94560', lw=1.2, alpha=0.85)
    ax1.axhline(0, color='#4A4A6A', lw=0.6, ls='--')
    ax1.axvline(0, color='#4A4A6A', lw=0.6, ls='--')
    ax1.set_xlabel("Campo B  [mT]",    color='#BDC3C7')
    ax1.set_ylabel("Magnetização  M (projeção)",  color='#BDC3C7')
    ax1.set_title("Curva de Histerese  M × B", color='#ECF0F1',
                  fontfamily='monospace')
    ax1.tick_params(colors='#7F8C8D')
    ax1.grid(True, color='#2C3E50', alpha=0.5)
    for sp in ax1.spines.values():
        sp.set_edgecolor('#2C3E50')

    # ── right panel: M(t) and B(t) versus time ──────────────────────────
    ax2 = axes[1]
    ax2.set_facecolor('#0D1B2A')
    l1, = ax2.plot(t_arr, M_arr, color='#E94560', lw=1.2, label='M (proj.)')
    ax2b = ax2.twinx()
    ax2b.plot(t_arr, B_arr, color='#FFD700', lw=1.0, alpha=0.7,
              ls='--', label='B (mT)')
    ax2b.set_ylabel("Campo B  [mT]", color='#FFD700')
    ax2b.tick_params(colors='#FFD700')
    ax2.set_xlabel("Tempo físico  t  [s]", color='#BDC3C7')
    ax2.set_ylabel("Magnetização  M", color='#E94560')
    ax2.set_title("Evolução temporal", color='#ECF0F1', fontfamily='monospace')
    ax2.tick_params(colors='#7F8C8D')
    ax2.grid(True, color='#2C3E50', alpha=0.4)
    for sp in ax2.spines.values():
        sp.set_edgecolor('#2C3E50')
    lines = [l1, plt.Line2D([0],[0], color='#FFD700', ls='--')]
    ax2.legend(lines, ['M (proj.)', 'B (mT)'], facecolor='#1A1A2E',
               edgecolor='#2C3E50', labelcolor='#ECF0F1', fontsize=8)

    fig.suptitle("Simulação de Histerese Magnética — Rede de Bússolas",
                 color='#BDC3C7', fontsize=12, fontfamily='monospace')
    plt.tight_layout()
    plt.savefig("hysteresis_loop.png", dpi=130, bbox_inches='tight',
                facecolor='#1A1A2E')
    plt.close(fig)
    _print("  Gráfico de histerese salvo: hysteresis_loop.png")


def _plot_sine(log, freq):
    """Plots M(t) and B(t) for the sinusoidal-field mode and saves the result as PNG.

Parameters
----------
log  : list of tuples (t, B_scalar, M_proj, S)
freq : sinusoidal-field frequency [Hz]"""
    log   = np.array(log)
    t_arr = log[:, 0]
    B_arr = log[:, 1] * 1e3    # T → mT
    M_arr = log[:, 2]

    fig, ax = plt.subplots(figsize=(10, 4), facecolor='#1A1A2E')
    ax.set_facecolor('#0D1B2A')
    ax.plot(t_arr, M_arr, color='#E94560', lw=1.2, label='M (proj.)')
    ax2 = ax.twinx()
    ax2.plot(t_arr, B_arr, color='#FFD700', lw=1.0, alpha=0.7,
             ls='--', label=f'B (mT)  f={freq:.2f} Hz')
    ax.set_xlabel("Tempo físico  t  [s]", color='#BDC3C7')
    ax.set_ylabel("Magnetização  M", color='#E94560')
    ax2.set_ylabel("Campo B  [mT]", color='#FFD700')
    ax.set_title(f"Campo Senoidal — f = {freq:.3f} Hz",
                 color='#ECF0F1', fontfamily='monospace')
    ax.tick_params(colors='#7F8C8D')
    ax2.tick_params(colors='#FFD700')
    ax.grid(True, color='#2C3E50', alpha=0.4)
    for sp in ax.spines.values():
        sp.set_edgecolor('#2C3E50')
    lines = [plt.Line2D([0],[0], color='#E94560'),
             plt.Line2D([0],[0], color='#FFD700', ls='--')]
    ax.legend(lines, ['M (proj.)', f'B (mT)  f={freq:.2f}Hz'],
              facecolor='#1A1A2E', edgecolor='#2C3E50',
              labelcolor='#ECF0F1', fontsize=8)
    plt.tight_layout()
    plt.savefig("sine_field.png", dpi=130, bbox_inches='tight',
                facecolor='#1A1A2E')
    plt.close(fig)
    _print("  Gráfico senoidal salvo: sine_field.png")


# ══════════════════════════════════════════════════════════════════════════════
# 3b. VECTORIZED CALCULATION OF THE DIPOLAR FIELD AND TORQUES (GPU/CPU)
# ══════════════════════════════════════════════════════════════════════════════

def compute_torques_vectorized(theta_flat, x_flat, y_flat, moment, cutoff,
                               bx_ext, by_ext,
                               pbc=False, Lx=None, Ly=None, n_images=1):
    """Computes the magnetic torque on each needle by summing the dipolar field from all other needles, fully vectorized without Python loops.

This function replaces the double loop `for i: for j: total_field_on(...)` with array operations (broadcasting), which makes it possible to run on the GPU through CuPy, whose syntax is identical to NumPy, when available, or normally on the CPU through NumPy. The complexity remains O(K²), where K = N×M, because every needle interacts with every other needle, but the sum is performed in parallel over all GPU/CPU cores at once instead of sequentially in Python. Typical speedups on GPUs such as the RTX 3090 are 10-100x for large lattices (K > about 400 needles).

Physics, identical to total_field_on/dipole_field_2d but vectorized:

    B = (mu0/4pi) · [ 3(m_hat·r_hat)r_hat - m_hat ] / r³

Parameters
----------
theta_flat : 1D array (K,) with all needle angles [rad]
x_flat, y_flat : 1D arrays (K,) with all needle positions [m]
moment     : magnetic moment of each needle [A·m²]
cutoff     : maximum interaction radius [m]
bx_ext, by_ext : instantaneous external-field components [T]
pbc        : if True, sums periodic replicas (see total_field_on)
Lx, Ly     : lattice periods in x and y [m], required if pbc=True
n_images   : number of replicas on each side in each direction (PBC)

Returns
-------
tau_flat : 1D array (K,) with the torque on each needle [N·m], in the same backend as theta_flat, GPU or CPU"""
    # detects the backend from the input array (instead of using the
    # module global _xp) — so the function correctly respects the
    # effective backend chosen by the caller (CPU or GPU via --gpu),
    # even if the global GPU detection is different
    if _GPU_AVAILABLE and hasattr(theta_flat, 'get'):
        xp = cp   # theta_flat is an array CuPy (GPU)
    else:
        xp = np   # theta_flat is an array NumPy (CPU)
    K  = theta_flat.shape[0]

    # magnetic-moment vectors of each source needle
    mx = moment * xp.cos(theta_flat)   # (K,)
    my = moment * xp.sin(theta_flat)   # (K,)

    # periodic shifts to sum (0.0 single one if pbc=False)
    if pbc and Lx and Ly:
        img_range = range(-n_images, n_images + 1)
        x_shifts = [k * Lx for k in img_range]
        y_shifts = [k * Ly for k in img_range]
    else:
        x_shifts = [0.0]
        y_shifts = [0.0]

    Bx_tot = xp.zeros(K)
    By_tot = xp.zeros(K)

    for dx_img in x_shifts:
        for dy_img in y_shifts:
            # distance matrix K×K: rx[the,b] = position(the) - position(b) - shift
            # broadcasting: x_flat[:, None] is a column (receiver), x_flat[None, :] is a row (source)
            rx = (x_flat[:, None] - x_flat[None, :]) - dx_img   # (K, K)
            ry = (y_flat[:, None] - y_flat[None, :]) - dy_img   # (K, K)

            r2 = rx*rx + ry*ry                                   # (K, K)

            # mask: distance > 0 (avoids self-interaction) and <= cutoff
            is_self = (r2 < 1e-24)
            valid   = (~is_self) & (r2 <= cutoff*cutoff)

            # avoids division by zero nos masked elements (will be zeroed afterwards)
            r2_safe = xp.where(valid, r2, 1.0)
            r       = xp.sqrt(r2_safe)
            r5      = r2_safe * r2_safe * r

            # m·r for each pair (the=receiver, b=source): uses mx[b], my[b] (source)
            mdotr = mx[None, :] * rx + my[None, :] * ry          # (K, K)

            bx_pair = MU0_OVER_4PI * (3.0 * mdotr * rx / r5 - mx[None, :] / (r2_safe * r))
            by_pair = MU0_OVER_4PI * (3.0 * mdotr * ry / r5 - my[None, :] / (r2_safe * r))

            # zeros invalid contributions (outside the cutoff or self-interaction)
            bx_pair = xp.where(valid, bx_pair, 0.0)
            by_pair = xp.where(valid, by_pair, 0.0)

            # sums over all sources (axis 1) for each receiver
            Bx_tot += bx_pair.sum(axis=1)
            By_tot += by_pair.sum(axis=1)

    Bx_tot += bx_ext
    By_tot += by_ext

    tau_flat = mx * By_tot - my * Bx_tot   # τ_z = mx·By − my·Bx
    return tau_flat


# ── V40: precomputed dipolar-interaction tensor ─────────────────────────
# Memory limit for the three blocks K×K of the tensor (float64). Above this,
# relax() automatically falls back for the per-step computestion (old method).
TENSOR_MEM_LIMIT_BYTES = 4.0e9   # 4 GB


def precompute_dipolar_tensor(x_flat, y_flat, cutoff,
                              pbc=False, Lx=None, Ly=None, n_images=1):
    """Precomputes the lattice dipolar-interaction tensor (V40).

Motivation: the positions of the needles are fixed throughout the simulation; only the angles theta evolve. In addition, the dipolar field is linear in the source moments (mx, my). Therefore, all pairwise geometry, including distances, powers of r, cutoff masks, and the sum over periodic images, can be condensed once before the integration loop into three constant K×K matrices:

    Bx_a = sum_b [ Axx[a,b]·mx_b + Axy[a,b]·my_b ]
    By_a = sum_b [ Axy[a,b]·mx_b + Ayy[a,b]·my_b ]   (Ayx = Axy)

with
    Axx[a,b] = (mu0/4pi)(3 rx²/r⁵ - 1/r³)
    Axy[a,b] = (mu0/4pi)(3 rx·ry/r⁵)
    Ayy[a,b] = (mu0/4pi)(3 ry²/r⁵ - 1/r³)

Each integration step is then reduced to matrix-vector products, BLAS on CPU or cuBLAS on GPU, instead of rebuilding about ten K×K arrays with square roots and masks at every step. The physics is the same, and pair values are the same up to floating-point reassociation; per-step cost is much lower.

The mask exactly reproduces compute_torques_vectorized semantics: self-interaction is excluded only in the zero-shift image, periodic replicas of the same needle contribute as they should, and the cutoff is applied image by image.

Parameters
----------
x_flat, y_flat : 1D arrays (K,) with positions [m], in the desired backend, NumPy/CPU or CuPy/GPU. The tensor is created in the same backend as the input arrays.
cutoff         : maximum interaction radius [m]
pbc, Lx, Ly, n_images : periodic-boundary settings, see relax()

Returns
-------
(Axx, Axy, Ayy) : three (K,K) matrices in the input backend, in units of T/(A·m²), field per unit moment."""
    if _GPU_AVAILABLE and hasattr(x_flat, 'get'):
        xp = cp
    else:
        xp = np
    K = x_flat.shape[0]

    if pbc and Lx and Ly:
        img_range = range(-n_images, n_images + 1)
        x_shifts = [k * Lx for k in img_range]
        y_shifts = [k * Ly for k in img_range]
    else:
        x_shifts = [0.0]
        y_shifts = [0.0]

    Axx = xp.zeros((K, K))
    Axy = xp.zeros((K, K))
    Ayy = xp.zeros((K, K))

    for dx_img in x_shifts:
        for dy_img in y_shifts:
            rx = (x_flat[:, None] - x_flat[None, :]) - dx_img
            ry = (y_flat[:, None] - y_flat[None, :]) - dy_img
            r2 = rx*rx + ry*ry

            is_self = (r2 < 1e-24)
            valid   = (~is_self) & (r2 <= cutoff*cutoff)

            r2_safe = xp.where(valid, r2, 1.0)
            r  = xp.sqrt(r2_safe)
            r3 = r2_safe * r
            r5 = r2_safe * r2_safe * r

            axx = MU0_OVER_4PI * (3.0 * rx * rx / r5 - 1.0 / r3)
            axy = MU0_OVER_4PI * (3.0 * rx * ry / r5)
            ayy = MU0_OVER_4PI * (3.0 * ry * ry / r5 - 1.0 / r3)

            Axx += xp.where(valid, axx, 0.0)
            Axy += xp.where(valid, axy, 0.0)
            Ayy += xp.where(valid, ayy, 0.0)

    return Axx, Axy, Ayy


def compute_torques_from_tensor(theta_flat, tensor, moment, bx_ext, by_ext):
    """Computes torques using the precomputed tensor (V40): four matrix-vector products per step. Physically identical to compute_torques_vectorized, with differences only from floating-point reassociation, around 1e-13 relative."""
    if _GPU_AVAILABLE and hasattr(theta_flat, 'get'):
        xp = cp
    else:
        xp = np
    Axx, Axy, Ayy = tensor
    mx = moment * xp.cos(theta_flat)
    my = moment * xp.sin(theta_flat)
    Bx = Axx @ mx + Axy @ my + bx_ext
    By = Axy @ mx + Ayy @ my + by_ext
    return mx * By - my * Bx


def relax(thetas, xs, ys, t_sim=2.0, dt_factor=0.05,
          inertia=INERTIA_DEFAULT, damping=DAMPING_DEFAULT,
          cutoff=3.5, ext_field=(0.0, 0.0), moment=MOMENT_DEFAULT,
          field_mode='static', field_freq=1.0,
          callback=None,
          frame_dir=None, frame_every=10,
          needle_len=0.042, needle_width=0.010, r_halo=None,
          frame_dpi=120, figsize_inches=None,
          pbc=False, Lx=None, Ly=None, n_images=1,
          B_ext=0.0, phi_ext_deg=0.0, use_gpu=False, show_progress=True,
          halo_mode='order', domain_tol_deg=15.0, make_images=True,
          hyst_spacing='linear', hyst_log_k=5.0,
          field_delay=0.0, t_pulse=None, torque_tol=1e-3,
          t_sim_full=False):
    """Integrates the real inertial dynamics of each needle with no pivot friction.

Physical model
--------------
Each needle is a rigid body with moment of inertia I rotating around the pivot, the z axis. The equation of motion is Newton's second law for rotation:

    I · d²theta/dt² = tau_z(theta, t) - b · dtheta/dt

where:
    I        [kg·m²] : moment of inertia, thin bar: I = mL²/12
    tau_z    [N·m]   : planar magnetic torque = mx·By - my·Bx
    b        [N·m·s] : viscous air damping
    dtheta/dt [rad/s]: angular velocity

Pivot friction = ZERO.

External-field modes (field_mode)
---------------------------------
'static'
    Constant field (Bext_x, Bext_y) over the whole interval.

'hysteresis'
    Field varies linearly in time over three ramps:
      0 -> t_sim/3          : B rises from 0 to +B_max
      t_sim/3 -> 2t_sim/3   : B falls from +B_max to -B_max
      2t_sim/3 -> t_sim     : B rises from -B_max to +B_max
    The direction is always phi_ext. The amplitude B_max is |ext_field|.
    The S=1.00 stop condition is disabled because magnetization oscillates.

'sine'
    Sinusoidal field: B(t) = B_max · sin(2*pi*f*t)
    Direction: phi_ext. Amplitude: B_max = |ext_field|. Frequency: field_freq [Hz].
    The S=1.00 stop condition is disabled.

Parameters
----------
thetas      : N×M array of initial angles [rad]
xs, ys      : fixed needle positions [m]
t_sim       : total physical simulation time [s], sum of dt steps, not CPU time
dt_factor   : fraction of the natural period T0 used as time step
inertia     : moment of inertia of each needle [kg·m²]
damping     : viscous damping coefficient b [N·m·s/rad]
cutoff      : maximum dipolar-interaction radius [m]
ext_field   : tuple (Bext_x, Bext_y), base field [T]. In 'static', constant field. In 'hysteresis'/'sine', defines direction and maximum amplitude.
moment      : magnetic moment of each needle [A·m²]
field_mode  : 'static' | 'hysteresis' | 'sine'
field_freq  : sinusoidal-field frequency [Hz], only for field_mode='sine'
callback    : optional function callback(step, thetas, omegas)
frame_dir   : directory for saving PNG frames, or None
frame_every : step interval between saved frames
make_images : bool, default True. If False, disables all image generation inside relax(): per-step PNG frames, even if frame_dir is passed, and final hysteresis_loop.png or sine_field.png summaries. Data CSVs and returned field_log are not affected.
needle_len  : needle size for rendering [m]
needle_width: needle width [m]
r_halo      : radius of order-parameter halos [m]
pbc         : if True, applies periodic boundary conditions in x and y so the lattice behaves as an infinite periodic structure
Lx, Ly      : lattice periods in x and y [m], required only if pbc=True and computed in make_grid
n_images    : number of periodic replicas to sum on each side in each direction when pbc=True; default 1 gives a (2·1+1)² = 9-cell replica grid
B_ext       : external-field magnitude for frame labels [T]
phi_ext_deg : external-field direction [degrees]

Returns
-------
theta_cur   : N×M array, final-state angles [rad]
omega_cur   : N×M array, angular velocities [rad/s]
hist        : list of tuples (thetas, omegas) every 20 steps
n_frames    : number of saved PNG frames
dt          : time step used [s]
stop_reason : string describing how the simulation ended"""
    import os
    import time as _time_module

    N, M = thetas.shape
    Bext_x, Bext_y = ext_field
    n_frames = 0

    # ── chooses the effective backend for this call (respects --gpu 0/1) ───────
    # _GPU_AVAILABLE (global of the module) indicates whether CuPy + GPU CUDA funcional
    # were detected at initialization. use_gpu (parameter, coming from --gpu)
    # decides whether this capability should actually be USED in this simulation.
    # This makes it possible to compare directly CPU vs GPU without restarting the process
    # or uninstalling the CuPy — useful for measuring the real performance gain.
    _active_gpu = use_gpu and _GPU_AVAILABLE
    _active_xp  = cp if _active_gpu else np

    def _local_to_backend(arr):
        """Local version of _to_backend that respects _active_gpu for this call."""
        if _active_gpu:
            return cp.asarray(arr)
        return np.asarray(arr)

    def _local_to_cpu(arr):
        """Local version of _to_cpu that respects _active_gpu for this call."""
        if _active_gpu and hasattr(arr, 'get'):
            return arr.get()
        return arr

    # ── nearest-neighbor distance ─────────────────────────────────
    # Vectorized via NumPy broadcasting (in time of 4 nested Python loops,
    # O(K²) interpreted iterations) — for large lattices (K = N×M in the
    # of the milhares), the pure-Python loop version podia levar many
    # segundos in this step alone, BEFORE the keyboard-listener thread (Ctrl+I/
    # Ctrl+C) is started — leaving the program "deaf" to any
    # interruption during that time. The vectorized version computes the same
    # coisa (minimum distance between any two distinct needles of the
    # lattice) in microseconds same for K~10000.
    x_flat_for_nn = xs.ravel()
    y_flat_for_nn = ys.ravel()
    dx_nn = x_flat_for_nn[:, None] - x_flat_for_nn[None, :]
    dy_nn = y_flat_for_nn[:, None] - y_flat_for_nn[None, :]
    d_nn  = np.sqrt(dx_nn**2 + dy_nn**2)
    np.fill_diagonal(d_nn, np.inf)   # excludes distance from each needle to itself
    r_nn  = float(np.min(d_nn))

    # ── effective field, natural frequency, dt and n_steps ───────────────────
    B_ref     = MU0_OVER_4PI * 2.0 * moment / r_nn**3
    B_ext_mag = np.sqrt(Bext_x**2 + Bext_y**2)
    B_eff     = max(B_ref, B_ext_mag)
    omega0    = np.sqrt(moment * B_eff / inertia)
    T0        = 2.0 * np.pi / omega0
    dt        = dt_factor * T0
    # number of steps computested from of the requested total time
    n_steps   = max(1, int(np.ceil(t_sim / dt)))
    Q         = omega0 * inertia / damping if damping > 0 else np.inf

    _print(f"  Dinamica inercial:")
    _print(f"    r_nn    = {r_nn*100:.2f} cm")
    _print(f"    B_ref   = {B_ref*1e3:.4f} mT  (dipolar entre vizinhos)")
    if B_ext_mag > 0:
        _print(f"    B_ext   = {B_ext_mag*1e3:.4f} mT  (campo externo)")
    _print(f"    B_eff   = {B_eff*1e3:.4f} mT  (campo dominante)")
    _print(f"    omega_0 = {omega0:.2f} rad/s   T0 = {T0:.5f} s")
    _print(f"    dt      = {dt:.6f} s  ({dt_factor:.0%} de T0)")
    _print(f"    t_sim   = {t_sim:.3f} s  -> {n_steps} passos")
    if Q > 2:
        q_desc = "sub-amortecido (oscila)"
    elif Q > 0.5:
        q_desc = "criticamente amortecido"
    else:
        q_desc = "super-amortecido"
    _print(f"    Q       = {Q:.1f}  ({q_desc})")
    _print()

    # ── initial conditions ─────────────────────────────────────────────────
    theta_cur = thetas.copy()
    omega_cur = np.zeros((N, M))
    hist      = [(theta_cur.copy(), omega_cur.copy())]

    # make_images=False turns off all image generation (PNG frames per step
    # and summary PNGs of hysteresis/sine), regardless of frame_dir ter
    # sido passado. Isso avoids matplotlib cost in programmatic sweeps
    # of parameter, without affecting the CSVs/field_log returned.
    if not make_images and frame_dir is not None:
        _print(f"  Aviso: make_images=False -> ignorando frame_dir='{frame_dir}' "
               f"(nenhum PNG sera gerado)")
        frame_dir = None

    if frame_dir is not None:
        os.makedirs(frame_dir, exist_ok=True)

    # ── formatting helper functions ──────────────────────────────────
    def _fmt_B(B):
        if B == 0:    return ""
        if B >= 0.1:  return f"B={B:.3f} T"
        if B >= 1e-4: return f"B={B*1e3:.3f} mT"
        return            f"B={B*1e6:.1f} µT"

    def _draw_clock(ax, t_phys, needle_len, stop_label=None):
        """Draws a stopwatch in the upper-left corner showing the system's physical time: the sum of integration steps dt, equivalent to the time a real stopwatch would measure while observing the needles move.

t_phys = step * dt [s]

Elements: main text with current time, total requested time, a progress bar, an early-stop label when applicable, and a semitransparent background box."""
        xlim  = ax.get_xlim()
        ylim  = ax.get_ylim()
        xspan = xlim[1] - xlim[0]
        yspan = ylim[1] - ylim[0]

        # panel position — upper-left corner
        px    = xlim[0] + 0.02 * xspan
        py    = ylim[1] - 0.03 * yspan

        bar_w = 0.30 * xspan
        bar_h = 0.018 * yspan
        bar_y = py - 0.095 * yspan

        extra_h = 0.030 * yspan if stop_label else 0.0

        # semi-transparent background box
        pad = needle_len * 0.25
        ax.add_patch(plt.Rectangle(
            (px - pad * 0.3, bar_y - pad * 0.6 - extra_h),
            bar_w + pad, 0.140 * yspan + pad + extra_h,
            facecolor='#080818', edgecolor='#3A3A6A',
            linewidth=0.8, alpha=0.82, zorder=19,
            transform=ax.transData))

        # ── progress bar: t_phys / t_sim ────────────────────────────
        frac = min(t_phys / t_sim, 1.0) if t_sim > 0 else 0.0

        # gray track
        ax.add_patch(plt.Rectangle(
            (px, bar_y), bar_w, bar_h,
            facecolor='#252540', edgecolor='none',
            zorder=20, transform=ax.transData))

        # fill: green → yellow → red
        if frac > 0:
            r_col = min(2.0 * frac, 1.0)
            g_col = min(2.0 * (1.0 - frac), 1.0)
            ax.add_patch(plt.Rectangle(
                (px, bar_y), bar_w * frac, bar_h,
                facecolor=(r_col, g_col, 0.15), edgecolor='none',
                zorder=21, transform=ax.transData))

        # border
        ax.add_patch(plt.Rectangle(
            (px, bar_y), bar_w, bar_h,
            facecolor='none', edgecolor='#5A5A8A',
            linewidth=0.7, zorder=22, transform=ax.transData))

        # white marker at the current position
        if 0 < frac < 1.0:
            mx = px + bar_w * frac
            ax.plot([mx, mx], [bar_y, bar_y + bar_h],
                    color='white', lw=1.2, zorder=24,
                    transform=ax.transData)

        # ── text: current physical time ──────────────────────────────────────
        # formats as mm:ss.ss if >= 60 s, otherwise in seconds with 4 decimal places
        if t_phys >= 60.0:
            mins = int(t_phys // 60)
            secs = t_phys - mins * 60
            time_str = f"t = {mins:02d}:{secs:05.2f}"
        else:
            time_str = f"t = {t_phys:.4f} s"

        ax.text(px + pad * 0.2, py,
                time_str,
                color='#E8E8FF', fontsize=11, fontweight='bold',
                fontfamily='monospace', va='top', ha='left',
                zorder=25, transform=ax.transData)

        # subtext: requested total time
        ax.text(px + pad * 0.2, py - 0.038 * yspan,
                f"/ {t_sim:.4f} s  (tempo físico)",
                color='#5555AA', fontsize=6,
                fontfamily='monospace', va='top', ha='left',
                zorder=25, transform=ax.transData)

        # stop label
        if stop_label:
            ax.text(px + pad * 0.2,
                    bar_y - pad * 0.5 - extra_h * 0.3,
                    stop_label,
                    color='#FFD700', fontsize=7, fontweight='bold',
                    fontfamily='monospace', va='top', ha='left',
                    zorder=25, transform=ax.transData)

    # ── field/magnetization history (initialized BEFORE of the functions of
    # frame, because the initial frame already queries hysteresis_log for the panel
    # M×H of the hysteresis) ──────────────────────────────────────────────────
    # field_log: (t, B, M_proj, S) for all modes.
    field_log      = []
    hysteresis_log = [] if field_mode == 'hysteresis' else None
    sine_log       = [] if field_mode == 'sine'       else None

    def _save_frame(step, th, om, stop_label=None,
                    B_ext_inst=None, phi_ext_inst=None):
        """Renders and saves a PNG frame.

B_ext_inst is the signed value of the field projected along phi [T]. It is negative during the reversed hysteresis phase. The panel shows the inverted arrow and signed value.

phi_ext_inst is the base field direction [degrees]. The closure variable B_ext is the maximum amplitude and is used to scale the bar.

V45: in hysteresis mode, creates a two-panel frame: needles on the left and the evolving M-H curve on the right as the field is swept."""
        nonlocal n_frames
        t_phys   = step * dt
        S        = np.abs(np.mean(np.exp(1j * th)))
        om_max   = np.max(np.abs(om))

        # signed field: B_now can be negative (reversed field)
        b_now    = B_ext_inst   if B_ext_inst   is not None else B_ext
        phi_now  = phi_ext_inst if phi_ext_inst is not None else phi_ext_deg

        # signed text: "+60.8 mT" or "-60.8 mT" or "0 T"
        if b_now is not None and abs(b_now) > 1e-12:
            sign_str = "+" if b_now > 0 else "-"
            b_str = sign_str + _fmt_B(abs(b_now))
        else:
            b_str = "0 T"

        title = f"S = {S:.4f}   w_max = {om_max:.2f} rad/s   B = {b_str}"

        # when B < 0: the arrow points in the opposite direction (phi + 180°)
        phi_display = phi_now if (b_now is None or b_now >= 0) else (phi_now + 180.0) % 360.0

        if field_mode == 'hysteresis':
            # ── two-panel frame: needles + curve M×H evolving ────────
            fig, (axL, axR) = plt.subplots(
                1, 2, figsize=(figsize_inches[0] * 2.05, figsize_inches[1]),
                facecolor='#1A1A2E',
                gridspec_kw={'width_ratios': [1.0, 1.0]})

            # left panel: needles (halos + compasss + field arrow)
            axL.set_facecolor('#16213E')
            axL.set_aspect('equal')
            _margin = needle_len * 1.6
            axL.set_xlim(xs.min() - _margin, xs.max() + _margin)
            axL.set_ylim(ys.min() - _margin, ys.max() + _margin * 2.5)
            axL.tick_params(left=False, bottom=False,
                            labelleft=False, labelbottom=False)
            for sp in axL.spines.values():
                sp.set_edgecolor('#2C3E50')
            draw_halos_batch(axL, xs, ys, th, needle_len, r_halo=r_halo,
                             halo_mode=halo_mode, domain_tol_deg=domain_tol_deg)
            draw_compass_batch(axL, xs, ys, th,
                               length=needle_len, width=needle_width)
            draw_ext_field_on_lattice(
                axL, xs, ys,
                abs(b_now) if b_now is not None else 0.0,
                phi_display, needle_len)
            axL.set_title(title, color='#ECF0F1', fontsize=11,
                          fontfamily='monospace')
            _draw_clock(axL, t_phys, needle_len, stop_label=stop_label)

            # right panel: curve M×H up to the current instant
            _render_hyst_panel(axR)

            plt.tight_layout()
            fpath = os.path.join(frame_dir, f"frame_{n_frames:05d}.png")
            plt.savefig(fpath, dpi=frame_dpi, bbox_inches='tight',
                        facecolor='#1A1A2E')
            plt.close(fig)
            n_frames += 1
            return

        fig, ax = plot_state(th, xs, ys, title=title,
                             needle_len=needle_len, needle_width=needle_width,
                             r_halo=r_halo,
                             B_ext=abs(b_now) if b_now is not None else 0.0,
                             phi_ext_deg=phi_display,
                             B_ext_max=B_ext,
                             B_signed=b_now,        # signed value for the text of the panel
                             figsize_inches=figsize_inches,
                             halo_mode=halo_mode,
                             domain_tol_deg=domain_tol_deg)
        plt.tight_layout()
        _draw_clock(ax, t_phys, needle_len, stop_label=stop_label)

        fpath = os.path.join(frame_dir, f"frame_{n_frames:05d}.png")
        plt.savefig(fpath, dpi=frame_dpi, bbox_inches='tight', facecolor='#1A1A2E')
        plt.close(fig)
        n_frames += 1

    def _render_hyst_panel(axR):
        """Draws the accumulated M-H curve up to the current step. Uses hysteresis_log, with the current point highlighted. Called at each frame in hysteresis mode (V45)."""
        axR.set_facecolor('#16213E')
        for sp in axR.spines.values():
            sp.set_edgecolor('#2C3E50')
        axR.tick_params(colors='#7F8C9A', labelsize=8)

        # fixed scale of H: ±B_max (in mT for readability)
        H_max_mT = (B_ext * 1e3) if B_ext > 0 else 1.0
        axR.set_xlim(-1.05 * H_max_mT, 1.05 * H_max_mT)
        axR.set_xlabel("H  (mT)", color='#BDC3C7', fontsize=10,
                       fontfamily='monospace')
        axR.set_ylabel("M  (projeção normalizada)", color='#BDC3C7',
                       fontsize=10, fontfamily='monospace')
        axR.set_title("Curva de histerese  M × H", color='#ECF0F1',
                      fontsize=11, fontfamily='monospace')

        # cross axes
        axR.axhline(0.0, color='#2C3E50', lw=1.0, zorder=1)
        axR.axvline(0.0, color='#2C3E50', lw=1.0, zorder=1)

        if hysteresis_log:
            arr = np.array(hysteresis_log)      # (t, B, M, S)
            H_mT = arr[:, 1] * 1e3
            Mv   = arr[:, 2]
            # M_proj is already the mean magnetization PER NEEDLE projetada na
            # field direction, normalized (saturates in ±1 when all
            # needles align) — does not require another normalization.
            axR.set_ylim(-1.15, 1.15)
            axR.plot(H_mT, Mv, color='#F1C40F', lw=1.6, zorder=3)
            # current point highlighted
            axR.plot([H_mT[-1]], [Mv[-1]], 'o', color='#E74C3C',
                     ms=7, zorder=4)

    if frame_dir is not None:
        _save_frame(0, theta_cur, omega_cur)

    # ── keyboard-listener thread (Ctrl+I = interactive interruption) ──────
    # Ctrl+I in the terminal sends the ASCII character \t (TAB, code 9).
    # A separate thread reads stdin in raw mode and sets the shared flag
    # when it detects \t, without interrupting the main process.
    # Works on Linux and macOS. On Windows it uses msvcrt.
    import threading, sys, os as _os

    _stop_flag = threading.Event()   # set by the thread or by S=1

    # original terminal state — saved in the outer scope for garantir
    # restoration even if the simulation ends for another reason other than
    # Ctrl+I (S=1, t_sim reached, rest, pulse completed). Without this, the
    # terminal may remain in raw mode after the program exits, deixando
    # the shell (zsh/bash) with broken output in subsequent commands.
    _term_fd  = None
    _term_old = None
    try:
        import termios
        _term_fd  = sys.stdin.fileno()
        _term_old = termios.tcgetattr(_term_fd)
    except Exception:
        pass   # stdin is not the terminal (pipe, redirection, etc.)

    def _restore_terminal():
        """Restores the terminal to normal cooked mode if it was changed."""
        if _term_fd is not None and _term_old is not None:
            try:
                import termios
                termios.tcsetattr(_term_fd, termios.TCSADRAIN, _term_old)
            except Exception:
                pass

    # registers as the final safety net: ensures restoration even in
    # case of an unhandled exception, Ctrl+C, or any abrupt termination
    import atexit
    atexit.register(_restore_terminal)

    def _keyboard_listener():
        """Keyboard-listener thread for Ctrl+I (Tab).

Goal: allow interrupting a long simulation and immediately finishing/saving the video, instead of losing all frames already generated.

Implementation:
- on POSIX systems, puts stdin in cbreak mode and reads one character at a time;
- Ctrl+I and Tab have ASCII code 9, so either key sets the stop flag;
- in non-interactive environments, the listener simply returns."""
        try:
            import tty, termios
            fd = sys.stdin.fileno()
            tty.setraw(fd)
            while not _stop_flag.is_set():
                ch = _os.read(fd, 1)
                if ch == b'\t':          # Ctrl+I = TAB = ASCII 9
                    _stop_flag.set()
                    break
                elif ch == b'\x03':      # Ctrl+C = ETX = ASCII 3
                    # IMMEDIATE and unconditional termination. In mode raw,
                    # the terminal does not automatically generate SIGINT from
                    # of Ctrl+C (this is disabled together with the rest of the
                    # special processing of control characters).
                    # Manually resending SIGINT (the.kill) was tested and
                    # proved unreliable: when the thread main
                    # is in a blocking call (e.g., waiting for this
                    # same thread, or inside a long NumPy operation),
                    # the signal delivery may never occur. Therefore here
                    # we terminate the process directly — without waiting the
                    # thread main "to notice" anything. Restauramos the
                    # terminal first, since the._exit() skips any
                    # finally/atexit from the rest of the program.
                    _restore_terminal()
                    sys.stdout.write("\r\n  Abortado (Ctrl+C) — encerramento imediato.\r\n")
                    sys.stdout.flush()
                    _os._exit(130)   # 130 = conventional exit code for SIGINT
        except Exception:
            # stdin is not the terminal (redirection, pipe, Windows without
            # msvcrt) — silently disables the keyboard listening
            pass
        finally:
            _restore_terminal()

    # starts thread the daemon so the not to block the program termination
    _kb_thread = threading.Thread(target=_keyboard_listener, daemon=True)
    _kb_thread.start()

    # ── base external field: direction and maximum amplitude ────────────────────
    Bext_x0, Bext_y0 = ext_field        # base field (direction + amplitude)
    B_max = np.sqrt(Bext_x0**2 + Bext_y0**2)   # maximum amplitude [T]
    phi_rad = np.arctan2(Bext_y0, Bext_x0)      # direction [rad]
    cos_phi = np.cos(phi_rad)
    sin_phi = np.sin(phi_rad)

    def field_at(t):
        """Computes the instantaneous signed external field for the selected field mode.

Returns (bx, by, B_scalar), where B_scalar is the signed projection along phi_ext. In static mode, the scalar is positive. In hysteresis and sine modes, it changes sign according to the ramp or sinusoid."""
        if field_mode == 'static':
            return Bext_x0, Bext_y0

        elif field_mode == 'hysteresis':
            # Full cycle: 0→+Bmax→0→-Bmax→0→+Bmax
            # 5 segments of t5 = t_sim/5. Within each segment, u0 is the
            # fraction traversed from of the endpoint B=0 of the segmento.
            #
            # hyst_spacing='linear' (default): |B| = Bmax·u0  (V40 and earlier)
            # hyst_spacing='log'   : |B| = Bmax·sinh(k·u0)/sinh(k)
            #     symmetric log-like spacing (via sinh, because the pure log
            #     does not cross zero): dB/dt small near B=0 → amostragem
            #     FINE field sampling in the transition/coercivity region; dB/dt
            #     large near ±Bmax → COARSE sampling at saturation.
            #     k (hyst_log_k) controls the concentratestion: k→0 recupera the
            #     linear; k=5 concentrates ~15x more points near of B=0 of the
            #     than near saturation.
            T  = t_sim if t_sim > 0 else 1.0
            t5 = T / 5.0
            if t <= t5:                      # seg 1: 0 → +Bmax
                u0, sgn = (t / t5), +1.0
            elif t <= 2.0 * t5:              # seg 2: +Bmax → 0
                u0, sgn = (1.0 - (t - t5) / t5), +1.0
            elif t <= 3.0 * t5:              # seg 3: 0 → -Bmax
                u0, sgn = ((t - 2.0 * t5) / t5), -1.0
            elif t <= 4.0 * t5:              # seg 4: -Bmax → 0
                u0, sgn = (1.0 - (t - 3.0 * t5) / t5), -1.0
            else:                            # seg 5: 0 → +Bmax
                u0, sgn = ((t - 4.0 * t5) / t5), +1.0
            if hyst_spacing == 'log' and hyst_log_k > 1e-9:
                # clamp only in log mode: the final step may fall beyond of
                # t_sim (u0 slightly > 1) and the sinh amplificaria esse
                # overshoot exponencialmente. No linear mode the overshoot
                # residual (~dt/t5 of B_max) is kept for compatibility
                # bit-for-bit with V40 and earlier.
                u0 = min(max(u0, 0.0), 1.0)
                g = np.sinh(hyst_log_k * u0) / np.sinh(hyst_log_k)
            else:
                g = u0
            B_scalar = sgn * B_max * g
            return B_scalar * cos_phi, B_scalar * sin_phi

        elif field_mode == 'sine':
            B_scalar = B_max * np.sin(2.0 * np.pi * field_freq * t)
            return B_scalar * cos_phi, B_scalar * sin_phi

        elif field_mode == 'pulse':
            # Timeline (V41):
            #   [0, field_delay)                : B = 0   (initial wait)
            #   [field_delay, end of the pulse)     : B = B_max
            #   [end of the pulse, ...)             : B = 0   (stabilization)
            # The "end of the pulse" is field_delay + t_pulse if t_pulse was given;
            # otherwise (legacy), it is when the criterion of S triggers
            # (_pulse_relaxing, controlled in the main loop).
            if t < field_delay:
                return 0.0, 0.0
            if t_pulse is not None:
                if t < field_delay + t_pulse:
                    return Bext_x0, Bext_y0
                return 0.0, 0.0
            if _pulse_relaxing[0]:
                return 0.0, 0.0
            return Bext_x0, Bext_y0

        elif field_mode == 'step_pos':
            # Positive step: B=0 during field_delay, then field
            # applied and held until the end.
            if t < field_delay:
                return 0.0, 0.0
            return Bext_x0, Bext_y0

        elif field_mode == 'step_neg':
            # Negative step: field applied from t=0; removed (zeroed)
            # in t=field_delay and remains so until the end.
            if t < field_delay:
                return Bext_x0, Bext_y0
            return 0.0, 0.0

        else:
            return Bext_x0, Bext_y0

    # mutable flag for the mode 'pulse' (one-element list for be
    # mutable inside _save_frame and the loop without nonlocal
    _pulse_relaxing = [False]   # False = field on; True = relaxing

    # prints usesge instruction
    _print("  [Ctrl+I para interromper e salvar o video agora]")
    _print("  [Ctrl+C para abortar IMEDIATAMENTE, sem salvar nada]")

    # ── helper function: computestes torques SI for state θ and field B(t) ────
    # ── pre-converts positions to "flat" arrays in the active backend ─────────
    # (GPU via CuPy if available, otherwise CPU via NumPy) — done only once
    # time before the integration loop, because xs/ys are fixed during all
    # the simulation (only the angles theta change at each step)
    x_flat_xp = _local_to_backend(xs.ravel())
    y_flat_xp = _local_to_backend(ys.ravel())

    # ── V40: precomputes the tensor of interaction dipolar ────────────────────
    # Fixed positions + field linear in the moments => all the geometry K×K
    # (distances, cutoff, images PBC) is cwherensed in 3 matrizes
    # constantes BEFORE the loop; each step becomes matrix-vector products
    # (BLAS/cuBLAS). Automatic fallback to the per-step method if the
    # tensor does not fit within the memory limit.
    _dipolar_tensor = None
    _K_total = int(x_flat_xp.shape[0])
    _tensor_bytes = 3.0 * _K_total * _K_total * 8.0
    if _tensor_bytes <= TENSOR_MEM_LIMIT_BYTES:
        _dipolar_tensor = precompute_dipolar_tensor(
            x_flat_xp, y_flat_xp, cutoff,
            pbc=pbc, Lx=Lx, Ly=Ly, n_images=n_images)
        _print(f"  Tensor dipolar pré-computado: 3 matrizes {_K_total}x{_K_total} "
               f"({_tensor_bytes/1e6:.0f} MB) — passos via matriz-vetor (V40)")
    else:
        _print(f"  Tensor dipolar NÃO pré-computado ({_tensor_bytes/1e9:.1f} GB "
               f"excede o limite de {TENSOR_MEM_LIMIT_BYTES/1e9:.0f} GB) — "
               f"usando cálculo por passo (método V39)")

    def _torques_xp(theta_flat_xp, bx_ext, by_ext):
        """GPU-resident version: receives and returns arrays in the active backend. This avoids CPU-GPU transfers inside the integration loop and is used when the simulation is actually running on the GPU."""
        if _dipolar_tensor is not None:
            return compute_torques_from_tensor(
                theta_flat_xp, _dipolar_tensor, moment, bx_ext, by_ext)
        return compute_torques_vectorized(
            theta_flat_xp, x_flat_xp, y_flat_xp, moment, cutoff,
            bx_ext, by_ext, pbc=pbc, Lx=Lx, Ly=Ly, n_images=n_images)

    def _torques(th, bx_ext, by_ext):
        """Compatibility version, CPU to GPU to CPU: accepts and returns CPU arrays, internally converting to the active backend only for the torque calculation. Useful for preserving older call sites and comparisons."""
        theta_flat_xp = _local_to_backend(th.ravel())
        tau_flat_xp = _torques_xp(theta_flat_xp, bx_ext, by_ext)
        tau_flat = _local_to_cpu(tau_flat_xp)
        return tau_flat.reshape(N, M)

    # ── performance instrumentation (wall-clock) ──────────────────────────
    # Measures real time spent in the integration loop (CPU or GPU), allowing
    # comparar diretamente the ganho of performance between --gpu 0 and --gpu 1.
    # Synchronizes the GPU before starting the timer, so that the measured time
    # reflect the actual work (without this, asynchronous GPU calls may
    # "voltar" before the computestion actually finishes, masking the real time).
    if _active_gpu:
        cp.cuda.Stream.null.synchronize()
    _perf_t_start  = _time_module.perf_counter()
    _perf_last_print_step = 0
    _perf_last_print_time = _perf_t_start
    _last_status_print_time = _perf_t_start
    _status_text = f"  [passo 0/{n_steps}]  aguardando primeiro status..."
    _last_bar_suffix = f"passo 0/{n_steps}"

    def _finish_progress_bar_if_shown():
        """Closes the progress bar only if it was enabled and visible. This avoids extra blank lines when progress output is disabled."""
        if show_progress:
            _print_progress_bar_finish()

    # ── Velocity-Verlet integration loop ────────────────────────────────
    # Main state (theta, omega, tau) resides in the active backend (_xp) —
    # GPU via CuPy if available, during the ENTIRE integration loop.
    # Isso avoids the overhead of transfer CPU<->GPU at each step, que
    # dominated runtime when only _torques() touched the GPU
    # (the rest of the loop, in NumPy puro, made constant round trips).
    # theta_cur/omega_cur (NumPy 2D) are kept the "mirrors" in CPU,
    # updated only when necessary (frame, history, callback).
    bx_cur, by_cur    = field_at(0.0)           # field at the initial instant
    theta_xp          = _local_to_backend(theta_cur.ravel())   # (K,) in the backend
    omega_xp           = _active_xp.zeros(N * M)                 # (K,) in the backend
    tau_xp             = _torques_xp(theta_xp, bx_cur, by_cur)
    _converged_count  = 0
    _S1_count         = 0   # consecutive steps with S >= 0.9999 (stop in static)
    # S=1 as a stopping criterion only makes sense in a static field.
    # In hysteresis and sine the field oscillates; in pulse we want to continue
    # after S=1 to observe the relaxation with zero field.
    # V42: t_sim_full=True disables ALL the early stops by
    # physical criterion (S=1.00, lattice in rest, equilibrium by torque) —
    # the simulation runs until the end of t_sim. User interruptions
    # (Ctrl+I, Ctrl+C) continue to work normally.
    _allow_S1_stop     = (field_mode == 'static') and (not t_sim_full)
    # in hysteresis and sine the field changes continuously: stop by rest
    # disabled (the lattice may be momentarily stopped at the B=0 crossing)
    # stop by rest: in 'static' always allowed; in 'pulse' only na
    # fase of relaxation (field already zeroed) — never during 'field_on', because
    # the lattice may remain momentarily stuck num frustrated local minimum
    # while the field is still on trying to align it
    _allow_rest_stop   = (field_mode == 'static') and (not t_sim_full)
    _stop_reason       = "tempo total atingido"
    # pulse phases (V41): 'delay' → 'field_on' → 'relaxing' → 'done'
    _pulse_phase   = "delay" if (field_mode == 'pulse' and field_delay > 0) \
                     else "field_on"
    _S99_count        = 0            # consecutive steps with S ≥ 0.99 (pulse)
    _S_window         = []           # sliding window of S recentes (pulse)
    # V41: reference torque for the equilibrium criterion of the modes
    # pulse/step_pos/step_neg: τ_ref = I·ω0² = m·B_eff (characteristic torque
    # of the dominant field). real equilibrium = ω_max small AND torque mean
    # |τ| < torque_tol·τ_ref (ω small alone may be only a point of
    # retorno of the oscillation, with large torque).
    _tau_ref = inertia * omega0 * omega0
    _step_after_announced = False    # message control of the modes step
    # window covers at least 2 natural oscillation periods (2 * 1/dt_factor
    # steps), ensuring that full cycles of oscillation are captured
    # before declaring stability — avoids False triggering in picos/vales
    # passageiros of the underdamped system
    _S_WINDOW_SIZE    = max(40, int(2.0 / dt_factor))

    for step in range(1, n_steps + 1):

        t_now = step * dt          # current physical time [s]

        # ── checks keyboard interruption (Ctrl+I) ──────────────────────
        if _stop_flag.is_set():
            _finish_progress_bar_if_shown()
            _print(f"\n  Interrompido (Ctrl+I) em t={t_now:.4f}s  (passo {step}/{n_steps})")
            _stop_reason = "interrompido pelo usuário (Ctrl+I)"
            if frame_dir is not None:
                theta_cur = _local_to_cpu(theta_xp).reshape(N, M)
                omega_cur = _local_to_cpu(omega_xp).reshape(N, M)
                _save_frame(step, theta_cur, omega_cur,
                            stop_label="■ interrompido (Ctrl+I)")
            break

        # ── external field at instant t (may vary with time) ──────────
        bx_new, by_new = field_at(t_now)

        # ── step 1: updates θ using θ(t), ω(t), τ(t) ───────────────────
        # (everything in "flat" arrays (K,) in the active backend — GPU if available)
        accel_xp  = (tau_xp - damping * omega_xp) / inertia
        theta_new_xp = theta_xp + omega_xp * dt + 0.5 * accel_xp * dt**2
        theta_new_xp = (theta_new_xp + _active_xp.pi) % (2.0 * _active_xp.pi) - _active_xp.pi

        # ── step 2: computes τ(t+dt) with the new θ and the novo field ──────────
        tau_new_xp = _torques_xp(theta_new_xp, bx_new, by_new)

        # ── step 3: updates ω (implicit Velocity-Verlet) ────────────────
        b_half    = damping * dt / (2.0 * inertia)
        omega_new_xp = (omega_xp * (1.0 - b_half)
                        + dt * (tau_xp + tau_new_xp) / (2.0 * inertia)) \
                       / (1.0 + b_half)

        theta_xp = theta_new_xp
        omega_xp = omega_new_xp
        tau_xp   = tau_new_xp
        bx_cur    = bx_new
        by_cur    = by_new

        # ── computes S, ω_max, mx_mean, my_mean ON THE GPU (reductions) ──────────
        # only the final scalars (4 floats) are brought to CPU the each
        # step, instead of the entire N×M arrays — drastically reducing the
        # volume of data transferred per integration step
        S_now_xp     = _active_xp.abs(_active_xp.mean(_active_xp.exp(1j * theta_xp)))
        omega_max_xp = _active_xp.max(_active_xp.abs(omega_xp))
        mx_mean_xp   = _active_xp.mean(_active_xp.cos(theta_xp))
        my_mean_xp   = _active_xp.mean(_active_xp.sin(theta_xp))

        S_now     = float(S_now_xp)
        omega_max = float(omega_max_xp)
        mx_mean   = float(mx_mean_xp)
        my_mean   = float(my_mean_xp)

        # ── syncs to 2D CPU only when necessary ────────────────
        # (callback, periodic history, frame): "lazy" conversion,
        # avoids the reshape/transfer cost when it is not needed
        _need_cpu_sync = (callback is not None) or (step % 20 == 0) or \
                         (frame_dir is not None and step % frame_every == 0)
        if _need_cpu_sync:
            theta_cur = _local_to_cpu(theta_xp).reshape(N, M)
            omega_cur = _local_to_cpu(omega_xp).reshape(N, M)

        if callback:
            callback(step, theta_cur.copy(), omega_cur.copy())

        if step % 20 == 0:
            hist.append((theta_cur.copy(), omega_cur.copy()))

        # ── records log of hysteresis / sine ──────────────────────────
        # M_proj: component of mean magnetization along the field direction
        # ── universal log: field and magnetization at each step ────────────────
        M_proj   = mx_mean * cos_phi + my_mean * sin_phi
        B_scalar = bx_cur * cos_phi + by_cur * sin_phi
        entry    = (t_now, B_scalar, M_proj, S_now)
        field_log.append(entry)
        if hysteresis_log is not None:
            hysteresis_log.append(entry)
        if sine_log is not None:
            sine_log.append(entry)

        # ── real-time progress bar (optional) ────────────────────
        # Visually shows the completed fraction of the integration, together with
        # throughput (steps/s) and the active backend (CPU/GPU). Updates at
        # maximum the each ~0.15s of wall-clock for the bar look smooth
        # without overloading the terminal with excessive writes. Pode be
        # disabled via --progress_bar 0 (ex: when redirecting output
        # for a log file, where the update through \r not se
        # comporta as in a interactive terminal).
        if show_progress:
            _perf_now = _time_module.perf_counter()
            if (_perf_now - _perf_last_print_time) >= 0.15 or step == n_steps:
                _steps_since = step - _perf_last_print_step
                _dt_wall     = _perf_now - _perf_last_print_time
                _steps_per_s = _steps_since / _dt_wall if _dt_wall > 0 else 0.0
                _backend_tag = "GPU" if _active_gpu else "CPU"
                _frac_done   = step / n_steps if n_steps > 0 else 1.0
                _last_bar_suffix = (f"passo {step}/{n_steps}  "
                                    f"({_steps_per_s:.0f} passos/s)")
                _print_progress_bar(
                    _frac_done,
                    prefix=f"  Integrando [{_backend_tag}] ",
                    suffix=_last_bar_suffix,
                    status_line=_status_text)
                _perf_last_print_step = step
                _perf_last_print_time = _perf_now

        # ── independent periodic status (field, S, ω) ────────────────────
        # Shows the evolution of the external field, order parameter and velocidade
        # angular over time, the each ~2s of wall-clock — ALWAYS active,
        # regardless of --progress_bar.
        # V44: with the bar enabled in an interactive terminal, the status is the
        # 2the row of the panel vivo, rewritten over itself (does not generate
        # a permanent line each cycle — which was what made the
        # bar "empilhar" rows). Without the bar, or in output redirecioanything
        # (not-TTY), keeps log behavior: a row nova by
        # ciclo via _print().
        _status_now = _time_module.perf_counter()
        if (_status_now - _last_status_print_time) >= 2.0 or step == n_steps:
            B_status = bx_cur * cos_phi + by_cur * sin_phi   # [T], with sign
            _status_text = (f"  [passo {step}/{n_steps}]  t={t_now:.4f}s  "
                            f"B={B_status*1e3:+.4f}mT  S={S_now:.4f}  "
                            f"w_max={omega_max:.3f}rad/s")
            if show_progress and _progress_ansi_ok():
                _backend_tag = "GPU" if _active_gpu else "CPU"
                _frac_done   = step / n_steps if n_steps > 0 else 1.0
                _print_progress_bar(
                    _frac_done,
                    prefix=f"  Integrando [{_backend_tag}] ",
                    suffix=_last_bar_suffix,
                    status_line=_status_text)
            else:
                if show_progress:
                    _print_progress_bar_finish()
                _print(_status_text)
            _last_status_print_time = _status_now

        # ── saves periodic frame ──────────────────────────────────────────
        if frame_dir is not None and step % frame_every == 0:
            # B_signed: signed projection along the field direction
            # (negative when the field points against phi_ext, ex: hysteresis)
            B_signed = bx_cur * cos_phi + by_cur * sin_phi   # [T], with sign
            _save_frame(step, theta_cur, omega_cur,
                        B_ext_inst=B_signed, phi_ext_inst=phi_ext_deg)
            # (progress now shown pela progress bar in time
            # real-time, above — a print per frame is no longer necessary)

        # ── condition 1: S = 1.00 (only in the static field) ──────────────────
        if _allow_S1_stop:
            if S_now >= 0.9999:
                _S1_count += 1
            else:
                _S1_count = 0

            if _S1_count >= 30:
                _finish_progress_bar_if_shown()
                _print(f"\n  S = 1.00 em t={t_now:.4f}s  (passo {step}/{n_steps})")
                _stop_reason = f"S = 1.00 atingido em t = {t_now:.4f} s"
                _stop_flag.set()
                if frame_dir is not None:
                    theta_cur = _local_to_cpu(theta_xp).reshape(N, M)
                    omega_cur = _local_to_cpu(omega_xp).reshape(N, M)
                    _save_frame(step, theta_cur, omega_cur,
                                stop_label="★ S = 1.00  alinhamento total")
                break

        # ── V41: timeral phase transitions (pulse with delay/t_pulse and steps) ──
        if field_mode == 'pulse':
            # end of the initial wait: field turns on
            if _pulse_phase == 'delay' and t_now >= field_delay:
                _pulse_phase = 'field_on'
                _finish_progress_bar_if_shown()
                _print(f"\n  Pulso: espera inicial concluída em t={t_now:.4f}s  campo LIGADO")
                tau_xp = _torques_xp(theta_xp, *field_at(t_now))
            # end of the pulse by fixed duration (t_pulse provided)
            if (_pulse_phase == 'field_on' and t_pulse is not None
                    and t_now >= field_delay + t_pulse):
                _pulse_relaxing[0] = True
                _pulse_phase = 'relaxing'
                _finish_progress_bar_if_shown()
                _print(f"\n  Pulso: duração t_pulse={t_pulse:.4f}s concluída em "
                       f"t={t_now:.4f}s  campo zerado, estabilizando")
                tau_xp = _torques_xp(theta_xp, 0.0, 0.0)
                if frame_dir is not None:
                    theta_cur = _local_to_cpu(theta_xp).reshape(N, M)
                    omega_cur = _local_to_cpu(omega_xp).reshape(N, M)
                    _save_frame(step, theta_cur, omega_cur,
                                stop_label="campo zerado - estabilizando",
                                B_ext_inst=0.0)
        elif field_mode in ('step_pos', 'step_neg') and not _step_after_announced:
            if t_now >= field_delay:
                _step_after_announced = True
                _finish_progress_bar_if_shown()
                if field_mode == 'step_pos':
                    _print(f"\n  Step+: campo LIGADO em t={t_now:.4f}s (permanece até o fim)")
                else:
                    _print(f"\n  Step-: campo REMOVIDO em t={t_now:.4f}s (permanece zerado)")
                tau_xp = _torques_xp(theta_xp, *field_at(t_now))

        # ── pulse mode (legacy, t_pulse=None): turns off the field when S>=0.99
        # OR when S stabilizes ──
        # Condition THE: S >= 0.99 for 20 consecutive steps (full alignment)
        # Condition B: inside a sliding window cobrindo pelo less
        #             2 natural oscillation periods (_S_WINDOW_SIZE
        #             steps), the AMPLITUDE of the S variation (max - min of the
        #             window) is smaller than 5% of the mean value. The wide window
        #             ensures that complete oscillation cycles are
        #             capturados, evitando disparo falso nos picos/vales
        #             passageiros of the underdamped system (where the
        #             instantaneous difference between steps may be small
        #             same longe of a stabilization real).
        if (field_mode == 'pulse' and _pulse_phase == 'field_on'
                and t_pulse is None):
            # condition THE: S >= 0.99
            if S_now >= 0.99:
                _S99_count += 1
            else:
                _S99_count = 0

            # condition B: amplitude of S inside of the sliding window < 5%
            _S_window.append(S_now)
            if len(_S_window) > _S_WINDOW_SIZE:
                _S_window.pop(0)

            _trigger_B = False
            if len(_S_window) >= _S_WINDOW_SIZE:
                S_max = max(_S_window)
                S_min = min(_S_window)
                S_mean = sum(_S_window) / len(_S_window)
                if S_mean > 1e-6:
                    spread_rel = (S_max - S_min) / S_mean
                    if spread_rel < 0.05:
                        _trigger_B = True

            _trigger_A = (_S99_count >= 20)

            if _trigger_A or _trigger_B:
                # arowmento reached OR S estabilizou → turns off the field
                _pulse_relaxing[0] = True
                _pulse_phase = 'relaxing'
                _finish_progress_bar_if_shown()
                if _trigger_A:
                    _print(f"\n  Pulso: S>=0.99 em t={t_now:.4f}s  campo zerado, relaxando")
                else:
                    _print(f"\n  Pulso: S estabilizou (S={S_now:.4f}) em t={t_now:.4f}s  campo zerado, relaxando")
                # updates tau with field=0 for the next step (keeps
                # the resident state is in the active backend — GPU if available)
                tau_xp = _torques_xp(theta_xp, 0.0, 0.0)
                if frame_dir is not None:
                    # B_ext_inst=0 shows the field off in the panel,
                    # but phi_ext_deg is kept to show the direction
                    theta_cur = _local_to_cpu(theta_xp).reshape(N, M)
                    omega_cur = _local_to_cpu(omega_xp).reshape(N, M)
                    _save_frame(step, theta_cur, omega_cur,
                                stop_label="campo zerado - relaxando",
                                B_ext_inst=0.0)

        # ── condition 2: lattice in equilibrium ─────────────────────────────────
        # In 'static': ω_max → 0 (comportamento original, inalterado).
        # In 'pulse' (fase 'relaxing') and nos modes 'step_pos'/'step_neg'
        # (after field_delay): V41 criterion for REAL equilibrium — ω_max
        # small AND torque mean |τ| < torque_tol·τ_ref, ambos by 50
        # steps. Torque is included because ω≈0 alone may be only a
        # turning point of the oscillation (with large torque); ω and τ are small
        # juntos caracterizam the equilibrium ("torque total small").
        # In 'hysteresis'/'sine': nunca, because the field changes continuously.
        _in_step_final = (field_mode in ('step_pos', 'step_neg')
                          and t_now >= field_delay)
        _rest_stop_active = (not t_sim_full) and (_allow_rest_stop or (
            field_mode == 'pulse' and _pulse_phase == 'relaxing') or
            _in_step_final)
        if _rest_stop_active:
            _needs_torque_check = (field_mode in ('pulse', 'step_pos',
                                                   'step_neg'))
            _eq_now = (omega_max < omega0 * 1e-3)
            if _eq_now and _needs_torque_check and torque_tol > 0:
                tau_mean = float(abs(tau_xp).mean())
                _eq_now = (tau_mean < torque_tol * _tau_ref)
            if _eq_now:
                _converged_count += 1
            else:
                _converged_count = 0

            if _converged_count >= 50:
                _finish_progress_bar_if_shown()
                _print(f"\n  Rede em repouso em t={t_now:.4f}s  (passo {step}/{n_steps})  S={S_now:.4f}")
                _stop_reason = f"rede em repouso em t = {t_now:.4f} s"
                _stop_flag.set()
                if frame_dir is not None:
                    theta_cur = _local_to_cpu(theta_xp).reshape(N, M)
                    omega_cur = _local_to_cpu(omega_xp).reshape(N, M)
                    _save_frame(step, theta_cur, omega_cur,
                                stop_label="● rede em repouso")
                break

    # ensures that the progress bar seja fechada even if the loop terminou
    # normalmente (atingiu n_steps) without passar by nenhuma of the mensagens
    # above, which already close the bar explicitly
    _finish_progress_bar_if_shown()

    # ── ensures final CPU<->GPU synchronization before returning ───────────
    # (the loop may have ended on a step where _need_cpu_sync was False)
    theta_cur = _local_to_cpu(theta_xp).reshape(N, M)
    omega_cur = _local_to_cpu(omega_xp).reshape(N, M)

    # ── final performance report (wall-clock) ─────────────────────────
    # Synchronizes the GPU before stopping the timer, so that the measured time
    # reflect the actual completed work (CUDA calls are asynchronous by
    # default — without this synchronization, the measured time could be smaller
    # than the real calculation time).
    if _active_gpu:
        cp.cuda.Stream.null.synchronize()
    _perf_t_end       = _time_module.perf_counter()
    _perf_total_s     = _perf_t_end - _perf_t_start
    _perf_steps_done  = step   # last step value executed in the loop
    _perf_steps_per_s = _perf_steps_done / _perf_total_s if _perf_total_s > 0 else 0.0
    _perf_ms_per_step = (_perf_total_s / _perf_steps_done * 1000.0) if _perf_steps_done > 0 else 0.0
    _perf_backend_tag = "GPU" if _active_gpu else "CPU"
    _print()
    _print(f"  ── desempenho ({_perf_backend_tag}) ──────────────────────────")
    _print(f"  Tempo de integracao : {_perf_total_s:.3f} s  (wall-clock)")
    _print(f"  Passos executados   : {_perf_steps_done}")
    _print(f"  Throughput          : {_perf_steps_per_s:.1f} passos/s  "
           f"({_perf_ms_per_step:.4f} ms/passo)")
    _print(f"  Agulhas na rede     : {N*M}  (K = N x M)")

    _stop_flag.set()   # ensures that the thread of keyboard encerra

    # restores the terminal explicitly here — not relying only on the thread,
    # because it may not have time to process _stop_flag before the
    # programa encerrar, deixando the terminal in raw mode for the shell
    _restore_terminal()

    # ── exporta dados of hysteresis / sine for CSV ─────────────────────
    if hysteresis_log:
        import csv
        csv_path = "hysteresis_loop.csv"
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['t_s', 'B_T', 'M_proj', 'S'])
            w.writerows(hysteresis_log)
        _print(f"  Dados de histerese salvos: {csv_path}")
        if make_images:
            _plot_hysteresis(hysteresis_log)

    if sine_log:
        import csv
        csv_path = "sine_field.csv"
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['t_s', 'B_T', 'M_proj', 'S'])
            w.writerows(sine_log)
        _print(f"  Dados de campo senoidal salvos: {csv_path}")
        if make_images:
            _plot_sine(sine_log, field_freq)

    return theta_cur, omega_cur, hist, n_frames, dt, _stop_reason, field_log

def next_available_path(path):
    """Builds a file path that does not overwrite an existing file by adding a numeric suffix when needed, for example name0001.mp4, name0002.mp4, and so on."""
    import os, re

    base, ext = os.path.splitext(path)

    # remove final numeric suffix (4 digits) from the base, if present
    # "sim0008" → "sim",  "dominios0003" → "dominios",  "sim" → "sim"
    prefix = re.sub(r'\d{4}$', '', base)

    # finds the next number in the sequence by checking only the video file
    n = 0
    while True:
        candidate = f"{prefix}{n:04d}{ext}"
        if not os.path.exists(candidate):
            return candidate
        n += 1


def render_video(frame_dir, output_path, fps=24, crf=20, use_gpu=False):
    """Aggregates numbered PNG frames from frame_dir into an MP4 video using ffmpeg.

Tries hardware-accelerated H.264 through NVENC first when available, then software H.264 through libx264, and finally MPEG-4 as a compatibility fallback. The explicit -f mp4 option avoids older ffmpeg builds failing to infer the output format from the file extension."""
    import subprocess, shutil, os

    if not shutil.which('ffmpeg'):
        _print("AVISO: ffmpeg não encontrado. Instale com:")
        _print("  Ubuntu/Debian : sudo apt install ffmpeg")
        _print("  macOS (brew)  : brew install ffmpeg")
        _print("  conda         : conda install -c conda-forge ffmpeg")
        _print("  Windows       : https://ffmpeg.org/download.html")
        return False

    input_pattern = os.path.join(frame_dir, "frame_%05d.png")
    vf = 'pad=ceil(iw/2)*2:ceil(ih/2)*2'

    strategies = []

    # NVENC is offered as the first option only if: (1) the user requested
    # explicitly --gpu 1 (same philosophy used for the physical calculation:
    # GPU is used only when requested), AND (2) an NVIDIA GPU was actually
    # detected on this machine (via CuPy). This avoids wasting time trying
    # NVENC on machines without an NVIDIA GPU, or when the user chose CPU.
    if use_gpu and _GPU_AVAILABLE:
        strategies.append((
            "H.264 via NVENC (GPU)",
            ['-c:v', 'h264_nvenc', '-preset', 'p4', '-cq', str(crf),
             '-pix_fmt', 'yuv420p']))

    strategies += [
        ("H.264 (libx264) — CPU",
         ['-c:v', 'libx264', '-preset', 'fast', '-crf', str(crf),
          '-pix_fmt', 'yuv420p']),
        ("MPEG-4 (mpeg4) — fallback para ffmpeg antigo",
         ['-c:v', 'mpeg4', '-q:v', '5',
          '-pix_fmt', 'yuv420p']),
        ("codec padrão do ffmpeg",
         ['-pix_fmt', 'yuv420p']),
    ]

    _print(f"\nMontando vídeo MP4: {output_path}")
    _print(f"  Fonte : {frame_dir}/frame_*.png   FPS={fps}   CRF={crf}")

    for desc, codec_args in strategies:
        cmd = [
            'ffmpeg', '-y',
            '-framerate', str(fps),
            '-i', input_pattern,
        ] + codec_args + [
            '-vf', vf,
            '-f', 'mp4',
            output_path,
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            size_mb = os.path.getsize(output_path) / 1024**2
            _print(f"  Codec : {desc}")
            _print(f"  Salvo : {output_path}  ({size_mb:.1f} MB)")
            return True
        else:
            last_err = result.stderr.strip().splitlines()
            short_err = last_err[-1] if last_err else "(sem mensagem)"
            _print(f"  [{desc}] falhou: {short_err}")

    _print(f"\nERRO: ffmpeg não conseguiu gerar o vídeo. Erro completo:")
    _print(result.stderr[-600:])
    return False
    """Aggregates numbered PNG frames from frame_dir into an MP4 video using ffmpeg.

Tries three strategies in order and stops at the first one that works:

1. H.264 (libx264) with explicit -f mp4, best quality and maximum compatibility with modern players.
2. MPEG-4 (mpeg4) with explicit -f mp4, fallback for older ffmpeg builds such as conda macOS builds without libx264.
3. Native ffmpeg codec, without -c:v, as a last resort.

The -f mp4 argument forces the output format explicitly and works around a bug in older ffmpeg versions that do not infer the format from the output-file extension.

Parameters
----------
frame_dir   : directory containing frame_00000.png, frame_00001.png, etc.
output_path : destination MP4 path
fps         : video frames per second

Returns
-------
True if video generation succeeded, False otherwise."""
    import subprocess, shutil, os

    if not shutil.which('ffmpeg'):
        _print("AVISO: ffmpeg não encontrado. Instale com:")
        _print("  Ubuntu/Debian : sudo apt install ffmpeg")
        _print("  macOS (brew)  : brew install ffmpeg")
        _print("  conda         : conda install -c conda-forge ffmpeg")
        _print("  Windows       : https://ffmpeg.org/download.html")
        return False

    input_pattern = os.path.join(frame_dir, "frame_%05d.png")

    # video filter: ensures even dimensions (required for H.264/MPEG-4)
    vf = 'pad=ceil(iw/2)*2:ceil(ih/2)*2'

    # list of strategies: (description, extra codec arguments)
    strategies = [
        ("H.264 (libx264)",
         ['-c:v', 'libx264', '-preset', 'fast', '-crf', str(crf),
          '-pix_fmt', 'yuv420p']),
        ("MPEG-4 (mpeg4) — fallback para ffmpeg antigo",
         ['-c:v', 'mpeg4', '-q:v', '5',
          '-pix_fmt', 'yuv420p']),
        ("codec padrão do ffmpeg",
         ['-pix_fmt', 'yuv420p']),
    ]

    _print(f"\nMontando vídeo MP4: {output_path}")
    _print(f"  Fonte : {frame_dir}/frame_*.png   FPS={fps}   CRF={crf}")

    for desc, codec_args in strategies:
        cmd = [
            'ffmpeg', '-y',
            '-framerate', str(fps),
            '-i', input_pattern,
        ] + codec_args + [
            '-vf', vf,
            '-f', 'mp4',          # force output format — fixes bug in older versions
            output_path,
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            size_mb = os.path.getsize(output_path) / 1024**2
            _print(f"  Codec : {desc}")
            _print(f"  Salvo : {output_path}  ({size_mb:.1f} MB)")
            return True
        else:
            # extracts only the last error line to avoid polluting the terminal
            last_err = result.stderr.strip().splitlines()
            short_err = last_err[-1] if last_err else "(sem mensagem)"
            _print(f"  [{desc}] falhou: {short_err}")

    # all strategies failed — shows the full error from the last attempt
    _print(f"\nERRO: ffmpeg não conseguiu gerar o vídeo. Erro completo:")
    _print(result.stderr[-600:])
    return False


# ══════════════════════════════════════════════════════════════════════════════
# 5. GENERATION OF THE NEEDLE GRID
# ══════════════════════════════════════════════════════════════════════════════

def make_grid(N=8, M=8, geometry='square', noise=1.5, R=0.025):
    """Creates the positions and initial angles of an N×M grid of needles.

The unifying parameter is R [m], the radius of the circle enclosing each needle. Positions are computed so adjacent circles touch tangentially, with center-to-center distance 2R, ensuring a visually correct lattice without circle overlap in any geometry.

Parameters
----------
N, M      : number of needle rows and columns
geometry  : lattice type
noise     : random-noise amplitude in the initial angles [rad]; 0 means all needles point to +x, and pi means fully random orientation
R         : radius of the circle enclosing each needle [m], default 0.025 m = 2.5 cm
seed      : random-number-generator seed

Returns
-------
thetas    : N×M array of initial angles [rad]
xs, ys    : N×M arrays of needle-center positions [m]
Lx, Ly    : lattice periods for periodic boundary conditions [m]"""
    s3 = np.sqrt(3.0)

    # ── grid square ─────────────────────────────────────────────────────
    if geometry == 'square':
        d = 2.0 * R          # distance between neighbors
        xs = np.zeros((N, M))
        ys = np.zeros((N, M))
        for i in range(N):
            for j in range(M):
                xs[i, j] = j * d
                ys[i, j] = i * d
        thetas = noise * np.random.randn(N, M)
        Lx, Ly = M * d, N * d   # period of the lattice for PBC [m]
        return xs, ys, thetas, d, Lx, Ly

    # ── grid triangular equilateral ───────────────────────────────────────
    elif geometry == 'triangular':
        d      = 2.0 * R          # distance between neighbors
        dx_col = d                # spacing between columns
        dy_row = R * s3           # spacing between rows = R·√3 = d·(√3/2)
        offset = R                # odd-row offset = R = d/2
        xs = np.zeros((N, M))
        ys = np.zeros((N, M))
        for i in range(N):
            x_off = offset * (i % 2)
            for j in range(M):
                xs[i, j] = j * dx_col + x_off
                ys[i, j] = i * dy_row
        thetas = noise * np.random.randn(N, M)
        Lx, Ly = M * dx_col, N * dy_row   # period of the lattice for PBC [m]
        return xs, ys, thetas, d, Lx, Ly

    # ── grid honeycomb (honeycomb) ──────────────────────────────────────────
    elif geometry == 'honeycomb':
        # ── Geometria of the image ──────────────────────────────────────────
        # Lines COMPLETAS alternando with rows of MEIA DENSITY,
        # creating visible hexagonal holes.
        #
        # Line pair  (row=0,2,4,...): M needles, step 2R, offset 0
        #   x = 0, 2R, 4R, ..., (M-1)*2R
        #
        # Line odd(row=1,3,5,...): ~M/2 needles, step 4R, offset R
        #   x = R, 5R, 9R, ...
        #
        # Δy between rows consecutivas = R*√3
        #
        # Each needle in the complete row touches:
        #   - 2 neighbors in the same row (at distance 2R)
        #   - 1 neighbor in the odd row above  (at distance 2R)  ← B(-R, +dy)
        #   - 1 neighbor in the odd row below (at distance 2R)
        #
        # Check: A(0,0) → B(-R, dy):
        #   d = √(R² + 3R²) = 2R  ✓
        #
        # Fill strategy:
        #   Generates a larger lattice (N+4 row pairs) and clips the W×H rectangle.
        # ─────────────────────────────────────────────────────────────────
        dy  = R * np.sqrt(3.0)   # spacing vertical between rows [m]
        d   = 2.0 * R            # distance between neighbors [m]

        # Target rectangle dimensions
        W = (M - 1) * 2.0 * R   # width: M needles in a complete row
        H = (N - 1) * dy         # altura:  N rows with spacing dy

        # Generates lattice larger (padding of 2 rows/columns in each side)
        N_rows = (N + 4) * 2     # * 2 porque alternamos pair and odd
        x_start = -2.0 * 2 * R   # starts two needles to the left
        y_start = -2.0 * dy       # starts two rows below

        xs_list, ys_list = [], []
        for row in range(N_rows):
            y = y_start + row * dy
            if row % 2 == 0:
                # complete row: x = x_start, x_start+2R, x_start+4R, ...
                x = x_start
                while x <= W + 2 * 2 * R:
                    xs_list.append(x)
                    ys_list.append(y)
                    x += 2.0 * R
            else:
                # row meia: x = x_start+R, x_start+5R, x_start+9R, ...
                # offset of R relative to the row pair, step 4R
                x = x_start + R
                while x <= W + 2 * 2 * R:
                    xs_list.append(x)
                    ys_list.append(y)
                    x += 4.0 * R

        all_x = np.array(xs_list)
        all_y = np.array(ys_list)

        # Clip: keeps points inside the target rectangle with margin R
        margin = R * 0.99
        mask = ((all_x >= -margin) & (all_x <= W + margin) &
                (all_y >= -margin) & (all_y <= H + margin))
        clipped_x = all_x[mask]
        clipped_y = all_y[mask]

        n_pts = len(clipped_x)
        xs     = clipped_x.reshape(n_pts, 1)
        ys     = clipped_y.reshape(n_pts, 1)
        thetas = noise * np.random.randn(n_pts, 1)
        Lx, Ly = W, H   # period aproximado of the lattice for PBC [m]
        return xs, ys, thetas, d, Lx, Ly

    else:
        raise ValueError(f"Geometria desconhecida: '{geometry}'. "
                         "Use 'square', 'triangular' ou 'honeycomb'.")


# ══════════════════════════════════════════════════════════════════════════════
# 5b. IDENTIFICATION OF MAGNETIC DOMAINS
# ══════════════════════════════════════════════════════════════════════════════

def label_magnetic_domains(thetas, tol_deg=15.0):
    """Identifies magnetic domains: connected regions of neighboring needles whose angular difference remains within a tolerance tol_deg.

Unlike the local order parameter, shown as green/red halos and measuring only the average alignment of each needle with its neighbors, this function groups the whole lattice into regions. Each region receives a distinct integer label, analogous to real magnetic domains: areas where magnetization points essentially in the same direction, separated by domain walls where orientation changes abruptly beyond the tolerance.

Algorithm
---------
Union-Find (disjoint-set) over the N×M grid. Direct neighbors (up/down/left/right, the same neighborhood used for the local order parameter) are merged when their angular difference modulo 2pi is <= tol_deg. The result is then relabeled compactly from 0 to n_domains-1.

Returns
-------
labels    : N×M integer array with domain labels
n_domains : number of domains"""
    N, M = thetas.shape
    K = N * M
    parent = np.arange(K)

    def find(x):
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    tol_rad = np.deg2rad(tol_deg)

    def _angdiff(a, b):
        # minimum angular difference on the circle, always in [0, π]
        d = a - b
        return np.abs(np.arctan2(np.sin(d), np.cos(d)))

    idx_grid = np.arange(K).reshape(N, M)

    # vertical pairs (i, i+1) whose angular difference is within the tolerance
    if N > 1:
        d_vert = _angdiff(thetas[:-1, :], thetas[1:, :])
        mask_vert = d_vert <= tol_rad
        pairs_a = [idx_grid[:-1, :][mask_vert]]
        pairs_b = [idx_grid[1:, :][mask_vert]]
    else:
        pairs_a, pairs_b = [], []

    # pairs horizontais (j, j+1)
    if M > 1:
        d_horiz = _angdiff(thetas[:, :-1], thetas[:, 1:])
        mask_horiz = d_horiz <= tol_rad
        pairs_a.append(idx_grid[:, :-1][mask_horiz])
        pairs_b.append(idx_grid[:, 1:][mask_horiz])

    if pairs_a:
        flat_a = np.concatenate(pairs_a)
        flat_b = np.concatenate(pairs_b)
        for a, b in zip(flat_a.tolist(), flat_b.tolist()):
            union(a, b)

    labels_raw = np.array([find(k) for k in range(K)])
    _, labels_remapped = np.unique(labels_raw, return_inverse=True)
    n_domains = int(labels_remapped.max()) + 1 if K > 0 else 0
    return labels_remapped.reshape(N, M), n_domains


def _domain_color_palette(n_domains):
    """Generates visually distinct and deterministic colors for n_domains domains, so the same domain label maps to the same color across frames and calls. Uses matplotlib's 'hsv' colormap sampled at regular hue intervals. Works well for a few domains and remains usable for many domains until the palette becomes visually dense, which is inherent to any palette with dozens of categories.

Parameters
----------
n_domains : number of domains to color

Returns
-------
array (n_domains, 4) with RGBA colors"""
    if n_domains <= 1:
        return np.array([[0.55, 0.55, 0.85, 1.0]])
    hues = np.linspace(0.0, 1.0, n_domains, endpoint=False)
    # small shift (golden ratio) prevents neighboring domains in
    # numeric label (which also tends to reflect spatial neighbors, since
    # the union-find rotula in order of varredura) caiam in matizes nearby
    hues = (hues * 0.618034 + 0.15) % 1.0
    return cm.hsv(hues)


# ══════════════════════════════════════════════════════════════════════════════
# 6. DESENHO OF A NEEDLE OF COMPASS
# ══════════════════════════════════════════════════════════════════════════════

def draw_compass(ax, x, y, theta, length=0.42, width=0.10,
                 color_n='#FFFFFF', color_s='#2E6DB4',
                 edge='#1a1a1a', zorder=4):
    """Draws a traditional compass needle as a two-color diamond.

Diamond geometry in local coordinates, with the long axis along +x:
    vertex 0 (+half, 0)      -> north tip, positive pole
    vertex 1 (0, +half_w)    -> upper width
    vertex 2 (-half, 0)      -> south tip, negative pole
    vertex 3 (0, -half_w)    -> lower width

The diamond is split into two halves at the center (x, y):
    north, white : triangle [0, 1, center, 3]
    south, blue  : triangle [2, 1, center, 3]

After construction in local coordinates, rotation by theta and translation to (x, y) are applied.

Parameters
----------
ax      : matplotlib axis on which to draw
x, y    : needle center or pivot position
theta   : needle angle [rad]
length  : needle length
width   : needle width"""
    half   = length / 2.0
    half_w = width  / 2.0

    # vertices in coordeanythings locais (not rotacionados)
    pts_local = np.array([
        [ half,     0.0    ],   # 0: tip north
        [ 0.0,      half_w ],   # 1: width upper
        [-half,     0.0    ],   # 2: tip south
        [ 0.0,     -half_w ],   # 3: width lower
    ])

    # matriz of rotation 2×2 pelo angle theta
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s],
                  [s,  c]])
    # applies rotation and translation
    pts = (R @ pts_local.T).T + np.array([x, y])

    # north half: polygon [north_tip, upper_width, center, lower_width]
    north = plt.Polygon(
        [pts[0], pts[1], [x, y], pts[3]],
        closed=True, facecolor=color_n,
        edgecolor=edge, linewidth=0.5, zorder=zorder)

    # south half: polygon [south_tip, upper_width, center, lower_width]
    south = plt.Polygon(
        [pts[2], pts[1], [x, y], pts[3]],
        closed=True, facecolor=color_s,
        edgecolor=edge, linewidth=0.5, zorder=zorder)

    ax.add_patch(south)
    ax.add_patch(north)

    # pino central (axis of rotation)
    ax.plot(x, y, 'o', ms=2.0, color='#555555',
            markeredgecolor='#222222', markeredgewidth=0.4,
            zorder=zorder + 1)


def draw_compass_batch(ax, xs, ys, thetas, length=0.42, width=0.10,
                       color_n='#FFFFFF', color_s='#2E6DB4',
                       edge='#1a1a1a', zorder=4):
    """Vectorized version of draw_compass(): draws all needles in the lattice at once using two PatchCollection objects for the north and south halves and a single ax.scatter() call for the central pivots, instead of one ax.add_patch()/ax.plot() call per needle.

Performance motivation
----------------------
Each individual ax.add_patch() triggers _update_patch_limits(), which traverses the path geometry, including Bezier curves, to update the axis bounding box. This fixed per-object cost becomes the dominant bottleneck when generating each frame with thousands of needles, measured around 1.5-1.85 ms per needle. PatchCollection receives all polygons at once and updates limits in batch, eliminating the repeated cost. The visual result is essentially identical, but rendering is much faster."""
    half   = length / 2.0
    half_w = width  / 2.0
    x_flat = np.asarray(xs).ravel()
    y_flat = np.asarray(ys).ravel()
    th_flat = np.asarray(thetas).ravel()
    K = x_flat.shape[0]

    # vertices locais (not rotacionados), same for all the needles
    pts_local = np.array([
        [ half,     0.0    ],   # 0: tip north
        [ 0.0,      half_w ],   # 1: width upper
        [-half,     0.0    ],   # 2: tip south
        [ 0.0,     -half_w ],   # 3: width lower
    ])   # shape (4, 2)

    # matrizes of rotation for all the K needles of a time: shape (K, 2, 2)
    c = np.cos(th_flat)
    s = np.sin(th_flat)
    # R[k] = [[c_k, -s_k], [s_k, c_k]]
    R = np.empty((K, 2, 2))
    R[:, 0, 0] = c
    R[:, 0, 1] = -s
    R[:, 1, 0] = s
    R[:, 1, 1] = c

    # applies the rotation of each needle to all four local vertices at once:
    # pts_rot[k, v, :] = R[k] @ pts_local[v]   →  shape (K, 4, 2)
    pts_rot = np.einsum('kab,vb->kva', R, pts_local)
    # translada each needle for sua position (x_k, y_k)
    pts_rot[:, :, 0] += x_flat[:, None]
    pts_rot[:, :, 1] += y_flat[:, None]

    # builds the north/south polygons of each needle from the vertices already
    # rotated and translated (same topology as the original version):
    #   north = [ponta_norte, largura_sup, center, largura_inf]
    #   south   = [ponta_sul,   largura_sup, center, largura_inf]
    centers = np.stack([x_flat, y_flat], axis=1)   # shape (K, 2)

    north_patches = []
    south_patches = []
    for k in range(K):
        p0, p1, p2, p3 = pts_rot[k]   # ponta_n, larg_sup, ponta_s, larg_inf
        ctr = centers[k]
        north_patches.append(mpatches.Polygon(
            [p0, p1, ctr, p3], closed=True))
        south_patches.append(mpatches.Polygon(
            [p2, p1, ctr, p3], closed=True))

    pc_north = PatchCollection(
        north_patches, facecolor=color_n, edgecolor=edge,
        linewidths=0.5, zorder=zorder)
    pc_south = PatchCollection(
        south_patches, facecolor=color_s, edgecolor=edge,
        linewidths=0.5, zorder=zorder)
    ax.add_collection(pc_south)
    ax.add_collection(pc_north)

    # central pins: a single scatter call for all needles, instead of
    # K individual calls to ax.plot()
    scatter_pins = ax.scatter(
        x_flat, y_flat, s=4.0, c='#555555',
        edgecolors='#222222', linewidths=0.4,
        zorder=zorder + 1)

    return pc_north, pc_south, scatter_pins


# ══════════════════════════════════════════════════════════════════════════════
# 7. LATTICE-STATE VISUALIZATION
# ══════════════════════════════════════════════════════════════════════════════

def _fmt_B(B):
    """Formats field intensity in a readable unit: T, mT, or microtesla."""
    if B == 0:      return "0 T"
    if B >= 0.1:    return f"{B:.4f} T"
    if B >= 1e-4:   return f"{B*1e3:.4f} mT"
    return              f"{B*1e6:.2f} µT"


def draw_ext_field_on_lattice(ax, xs, ys, B_ext, phi_ext_deg,
                              needle_len, B_ext_max=None, B_signed=None,
                              color='#FFD700'):
    """Represents the external field in the image. It is always drawn, even when B_ext = 0, for example after a pulse has been turned off.

Elements
--------
1. ARROWS AT SITES: golden arrows at each needle pointing in the field direction. When B_ext = 0, the arrows become invisible (alpha=0) but the reference panel remains visible.

2. FIELD PANEL: box in the upper-right corner, mirroring the stopwatch geometry in the upper-left corner and avoiding overlap with the needles. It contains the current SI intensity and angle in degrees, a horizontal intensity bar, and a small direction arrow to the right of the bar. When B_ext = 0, the arrow is hidden, the text reads "0 T", and the bar is empty.

Parameters
----------
ax           : matplotlib axis
xs, ys       : needle positions
B_ext        : external-field magnitude [T]
phi_ext_deg  : external-field direction [degrees]
needle_len   : reference needle length for scaling"""
    phi      = np.deg2rad(phi_ext_deg)
    cos_phi  = np.cos(phi)
    sin_phi  = np.sin(phi)
    B_max_ref = B_ext_max if (B_ext_max and B_ext_max > 0) else max(B_ext, 1e-30)
    frac      = min(abs(B_ext) / B_max_ref, 1.0)   # 0→1, proportion of the maximum field
    field_on  = (B_ext > 1e-30)                     # True if the field is nonzero

    # ── 2. panel field (ALWAYS visible) — corner upper right ──────
    # Mirrors the clock geometry (upper-left corner): box
    # retangular ancorada ao top, text and bar horizontal insidela.
    # This avoids overlap with the needles, which remain in the central area.
    xlim  = ax.get_xlim()
    ylim  = ax.get_ylim()
    xspan = xlim[1] - xlim[0]
    yspan = ylim[1] - ylim[0]

    bar_w = 0.20 * xspan
    bar_h = 0.018 * yspan
    arrow_zone_w = 0.07 * xspan   # space reserved for the direction arrow

    # upper-right corner: anchored to the top, with a right margin
    box_w = bar_w + arrow_zone_w
    px = xlim[1] - 0.02 * xspan - box_w    # border left of the caixa
    py = ylim[1] - 0.03 * yspan             # top of the box (same level as the clock)

    pad = needle_len * 0.25

    # semi-transparent background box — same visual style as the clock
    box_h = 0.140 * yspan + pad
    ax.add_patch(plt.Rectangle(
        (px - pad * 0.3, py - box_h),
        box_w + pad, box_h,
        facecolor='#080818', edgecolor=color,
        linewidth=0.8, alpha=0.82, zorder=19,
        transform=ax.transData))

    # ── text: intensidade and direction (row upper of the caixa) ────────────
    b_val = B_signed if B_signed is not None else (B_ext if field_on else 0.0)
    if abs(b_val) > 1e-12:
        sign_str = "+" if b_val >= 0 else "-"
        b_str = sign_str + _fmt_B(abs(b_val))
    else:
        b_str = "0 T  (desligado)" if not field_on else "0 T"
    phi_str = f"dir: {phi_ext_deg:.1f} graus"

    text_color = color if field_on else '#777777'
    ax.text(px, py - 0.022 * yspan, "B ext", color=color, fontsize=7,
            alpha=0.8, fontfamily='monospace', ha='left', va='top', zorder=23)
    ax.text(px, py - 0.048 * yspan, b_str, color=text_color,
            fontsize=9, fontweight='bold', fontfamily='monospace',
            ha='left', va='top', zorder=23)
    ax.text(px, py - 0.072 * yspan, phi_str, color=text_color,
            fontsize=7, fontfamily='monospace', ha='left', va='top', zorder=23)

    # ── horizontal intensity bar (same style as the clock) ───
    bar_y = py - 0.115 * yspan
    # trilha
    ax.add_patch(plt.Rectangle(
        (px, bar_y), bar_w, bar_h,
        facecolor='#252540', edgecolor='none',
        zorder=20, transform=ax.transData))
    # fill proporcional ao field current / maximum
    if frac > 0.01:
        ax.add_patch(plt.Rectangle(
            (px, bar_y), bar_w * frac, bar_h,
            facecolor=color if field_on else '#555555',
            edgecolor='none', alpha=0.85, zorder=21,
            transform=ax.transData))
    # border
    ax.add_patch(plt.Rectangle(
        (px, bar_y), bar_w, bar_h,
        facecolor='none', edgecolor='#5A5A8A',
        linewidth=0.7, zorder=22, transform=ax.transData))

    # ── small direction arrow, to the right of the bar, inside the box ──────
    arrow_cx = px + bar_w + arrow_zone_w / 2.0
    arrow_cy = bar_y + bar_h / 2.0
    alen     = min(needle_len * 0.7, arrow_zone_w * 0.8)
    arrow_len = alen * max(frac, 0.0)
    if arrow_len > alen * 0.05:
        x_tail = arrow_cx - arrow_len * cos_phi / 2.0
        y_tail = arrow_cy - arrow_len * sin_phi / 2.0
        x_head = arrow_cx + arrow_len * cos_phi / 2.0
        y_head = arrow_cy + arrow_len * sin_phi / 2.0
        arrow_color = color if field_on else '#555555'
        ax.annotate(
            '', xy=(x_head, y_head), xytext=(x_tail, y_tail),
            arrowprops=dict(
                arrowstyle='->', color=arrow_color, lw=1.8,
                mutation_scale=12),
            zorder=20)
    else:
        ax.plot(arrow_cx, arrow_cy, 'o', ms=3, color='#555555', zorder=20)

    # ── 1. arrows at the sites (only when the field is on) ─────────────────
    # ax.quiver draws ALL arrows in a single call (instead of
    # N×M individual calls to ax.annotate, each creating an arrow patch
    # arrow patch with its own bounding-box update cost).
    # Because the field is uniform, all arrows have the same direction/magnitude
    # — quiver aceita arrays of position with escalares of component (U, V),
    # broadcasting automatically to all points.
    if field_on:
        alen = needle_len * 0.35
        dx_a = alen * cos_phi
        dy_a = alen * sin_phi
        N, M = xs.shape
        x0 = xs.ravel() - dx_a / 2.0
        y0 = ys.ravel() - dy_a / 2.0
        ax.quiver(x0, y0, dx_a, dy_a,
                 angles='xy', scale_units='xy', scale=1.0,
                 color=color, alpha=0.50, width=0.0035,
                 headwidth=4.0, headlength=4.5, headaxislength=4.0,
                 zorder=2)


def draw_halos_batch(ax, xs, ys, thetas, needle_len, r_halo=None,
                     halo_mode='order', domain_tol_deg=15.0):
    """Draws the colored halos of all needles on an axis (V43). The logic was extracted from plot_state for reuse in the side-by-side comparison figure and in other panels.

The calculation is vectorized through array shifts instead of Python loops over needles, and drawing uses PatchCollection, one call for all circles instead of N×M individual calls.

Modes (halo_mode):
  'order'   : green means aligned with neighbors and red means frustrated, based on the average cos(delta theta) with direct grid neighbors.
  'domains' : each connected region with similar orientation, within domain_tol_deg, receives a distinct categorical color; small domains, smaller than max(3, 2% of needles), are shown in neutral gray to reduce visual noise."""
    N, M = thetas.shape
    _r = r_halo if r_halo is not None else needle_len * 0.58
    _n_domains_found = None
    _n_significant   = None

    if halo_mode == 'domains':
        # ── domain mode: each connected region with similar orientation
        # receives a distinct categorical color, identifying magnetic domains
        # (domain walls = abrupt color change).
        #
        # Very small domains (typically one to a few needles, caused
        # from residual noise in a lattice that is not fully relaxed yet)
        # receive a color CINZA NEUTRA in time of a color distinta of the
        # palette — otherwise, with hundreds of small domains,
        # the categorical palette visually degenerates into a gradient
        # continuous (adjacent hues very close), hiding the
        # few large domains that actually matter for interpretation.
        # The threshold is proportional to the lattice size (or an absolute minimum
        # absolute minimum of 3 needles, whichever is larger), and the domains
        # "significactives" remanescentes receive colors ordeanythings of the
        # larger for the smaller, maximizando contraste between domains
        # spatially neighboring domains (which tend to have numerical labels
        # nearby, since the union-find rotula in order of varredura).
        domain_labels, _n_domains_found = label_magnetic_domains(
            thetas, tol_deg=domain_tol_deg)
        K_total = N * M
        _unique_labels, _inverse, _counts = np.unique(
            domain_labels, return_inverse=True, return_counts=True)
        _inverse = _inverse.ravel()   # Recent NumPy returns the same shape as the input (N,M); we need (K,)
        _min_domain_size = max(3, int(0.02 * K_total))
        _is_significant = _counts >= _min_domain_size
        _n_significant = int(np.sum(_is_significant))

        # sorts domains significactives of the larger for the smaller, for
        # so that the palette assigns colors in the same order (largest domain
        # receives the primeira color of the palette, and assim by diante)
        _sig_idx = np.where(_is_significant)[0]
        _sig_order = _sig_idx[np.argsort(-_counts[_sig_idx])]
        _color_rank = np.full(len(_unique_labels), -1, dtype=int)
        _color_rank[_sig_order] = np.arange(len(_sig_order))

        palette = _domain_color_palette(max(_n_significant, 1))
        _gray = np.array([0.45, 0.45, 0.48, 1.0])
        per_label_colors = np.tile(_gray, (len(_unique_labels), 1))
        if _n_significant > 0:
            per_label_colors[_sig_order] = palette[:_n_significant]

        colors = per_label_colors[_inverse]
    else:
        # ── order mode (default): green=aligned, red=frustrated —
        # mean of cos(Δθ) with the neighbors diretos (until 4: cima/baixo/
        # left/right), tratando borders of the grid corretamente
        # (neighbors outside of the grid not contam)
        align_sum   = np.zeros((N, M))
        align_count = np.zeros((N, M))
        for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            # dst = region of (i,j) that has a valid neighbor at (i+di, j+dj);
            # src = same region, shifted by (di, dj) — it is the corresponding neighbor
            # correspwherente of each elemento of dst, pair the pair
            dst_i0, dst_i1 = max(0, -di), N - max(0, di)
            dst_j0, dst_j1 = max(0, -dj), M - max(0, dj)
            src_i0, src_i1 = dst_i0 + di, dst_i1 + di
            src_j0, src_j1 = dst_j0 + dj, dst_j1 + dj
            delta = thetas[dst_i0:dst_i1, dst_j0:dst_j1] - thetas[src_i0:src_i1, src_j0:src_j1]
            align_sum[dst_i0:dst_i1, dst_j0:dst_j1]   += np.cos(delta)
            align_count[dst_i0:dst_i1, dst_j0:dst_j1] += 1.0

        align_avg = np.where(align_count > 0, align_sum / np.maximum(align_count, 1.0), 0.0)
        colors = cm.RdYlGn((align_avg.ravel() + 1.0) / 2.0)

    halo_patches = [plt.Circle((x, y), _r)
                    for x, y in zip(xs.ravel(), ys.ravel())]
    _halo_alpha = 0.35 if halo_mode == 'domains' else 0.20
    pc_halos = PatchCollection(
        halo_patches, facecolor=colors, edgecolor='none',
        alpha=_halo_alpha, zorder=1)
    ax.add_collection(pc_halos)
    return _n_domains_found, _n_significant


def plot_state(thetas, xs, ys, title="Rede de bussolas", show_order=True,
               needle_len=0.42, needle_width=0.10, r_halo=None,
               B_ext=0.0, phi_ext_deg=0.0, B_ext_max=None,
               B_signed=None, figsize_inches=None,
               halo_mode='order', domain_tol_deg=15.0):
    """Generates a figure with the instantaneous state of the needle lattice.

Visual elements
---------------
- Colored halo around each needle. Two modes are available:
    'order', default: green = neighbors aligned, red = anti-aligned/frustrated; local order parameter from the average cos(delta theta) with grid neighbors.
    'domains': each connected region of similarly oriented needles, within domain_tol_deg, receives a distinct categorical color, identifying magnetic domains.
- Needles drawn as two-color diamonds, white = north and blue = south.
- Golden arrows at each site showing the external-field direction.
- A lower-right panel with the field magnitude/direction.

Parameters
----------
thetas      : N×M array of needle angles
xs, ys      : N×M arrays of positions
outpath     : PNG output path
needle_len  : rendered needle length
needle_width: rendered needle width
B_ext       : external-field magnitude [T]
phi_ext_deg : external-field direction [degrees]
halo_mode   : 'order' or 'domains'
domain_tol_deg : angular tolerance for domain grouping [degrees]"""
    # ── size of the figure: escala with the extension physics of the lattice ──────────
    # Ensures needles have adequate visual size regardless of
    # how many there are in the lattice. The figure grows with the lattice but is limited.
    if figsize_inches is None:
        x_span = xs.max() - xs.min() + 2 * needle_len * 2.0
        y_span = ys.max() - ys.min() + 2 * needle_len * 4.0
        aspect = x_span / y_span if y_span > 0 else 1.0
        # base: 8 inches on the smaller side, maximum 20 inches
        base = 8.0
        if aspect >= 1.0:
            fig_w = min(base * aspect, 20.0)
            fig_h = fig_w / aspect
        else:
            fig_h = min(base / aspect, 20.0)
            fig_w = fig_h * aspect
        figsize_inches = (max(fig_w, 6.0), max(fig_h, 6.0))

    fig, ax = plt.subplots(figsize=figsize_inches, facecolor='#1A1A2E')
    ax.set_facecolor('#16213E')
    N, M = thetas.shape

    # ── halos coloridos: order parameter local OR domains magnetic ──
    # (logic extracted to draw_halos_batch — V43 — for reuse in the
    # side-by-side comparison figure; behavior identical to the previous one)
    if show_order:
        _n_domains_found, _n_significant = draw_halos_batch(
            ax, xs, ys, thetas, needle_len, r_halo=r_halo,
            halo_mode=halo_mode, domain_tol_deg=domain_tol_deg)

    # ── needles ───────────────────────────────────────────────────────────
    draw_compass_batch(ax, xs, ys, thetas,
                       length=needle_len, width=needle_width)

    # ── axis formatting (must come BEFORE the field panel) ──────────
    ax.set_aspect('equal')
    margin = needle_len * 1.6
    # larger upper margin: reserves space for the clock panels
    # (corner upper left) and external field (corner upper right), with
    # additional empty space below the panels until the first row
    # of needles, evitando any sobreposition
    top_margin = needle_len * 7.0
    ax.set_xlim(xs.min() - margin, xs.max() + margin)
    ax.set_ylim(ys.min() - margin, ys.max() + top_margin)
    ax.set_title(title, color='#ECF0F1', fontsize=11, pad=10,
                 fontfamily='monospace')
    ax.tick_params(colors='#7F8C8D')
    for spine in ax.spines.values():
        spine.set_edgecolor('#2C3E50')

    # ── external field (arrows + panel) — debecause of set_xlim/ylim ─────────
    draw_ext_field_on_lattice(ax, xs, ys, B_ext, phi_ext_deg, needle_len,
                              B_ext_max=B_ext_max, B_signed=B_signed)

    # ── legenda ───────────────────────────────────────────────────────────
    patch_n = mpatches.Patch(facecolor='#FFFFFF', edgecolor='#555',
                              label='Polo Norte (branco)')
    patch_s = mpatches.Patch(facecolor='#2E6DB4', edgecolor='#555',
                              label='Polo Sul (azul)')
    handles = [patch_n, patch_s]
    if show_order and halo_mode == 'domains':
        if _n_significant < _n_domains_found:
            _dom_label = (f'Dominios: {_n_significant} principais '
                         f'(+{_n_domains_found - _n_significant} pequenos, cinza)  '
                         f'tol={domain_tol_deg:.0f}°')
        else:
            _dom_label = f'Dominios magneticos: {_n_domains_found}  (tol={domain_tol_deg:.0f}°)'
        patch_dom = mpatches.Patch(facecolor='none', edgecolor='none',
                                   label=_dom_label)
        handles.append(patch_dom)
    if B_ext > 0:
        patch_b = mpatches.Patch(facecolor='#FFD700', edgecolor='#555',
                                  label=f'Campo ext.  {_fmt_B(B_ext)}  '
                                        f'φ={phi_ext_deg:.1f}°')
        handles.append(patch_b)
    ax.legend(handles=handles, loc='lower left',
              facecolor='#1A1A2E', edgecolor='#2C3E50',
              labelcolor='#ECF0F1', fontsize=8)

    return fig, ax


# ══════════════════════════════════════════════════════════════════════════════
# 8. RELAXATION ANIMATION
# ══════════════════════════════════════════════════════════════════════════════

def animate_relaxation(thetas_hist, xs, ys,
                       needle_len=0.42, needle_width=0.10,
                       B_ext=0.0, phi_ext_deg=0.0,
                       interval=80, save_gif=None):
    """Generates a frame-by-frame animation of the lattice time evolution.

Each frame corresponds to a snapshot saved during relaxation, one every 20 integration steps.

Parameters
----------
thetas_hist  : list of N×M arrays from the `hist` field returned by `relax`
xs, ys       : needle positions
needle_len   : needle size
needle_width : needle width
B_ext        : external-field magnitude, for the arrow
phi_ext_deg  : external-field direction, in degrees
interval     : time between frames in ms, smaller means faster
save_gif     : path for saving the GIF; None displays on screen

Returns
-------
fig, ani : matplotlib objects (FuncAnimation)"""
    fig, ax = plt.subplots(figsize=(7, 7), facecolor='#1A1A2E')
    ax.set_facecolor('#16213E')
    ax.set_aspect('equal')
    margin = needle_len * 1.6
    ax.set_xlim(xs.min() - margin, xs.max() + margin)
    ax.set_ylim(ys.min() - margin, ys.max() + margin * 2.5)

    N, M = thetas_hist[0].shape

    def update(frame):
        # clears all artists from the previous frame
        for art in ax.lines + ax.patches + ax.collections:
            art.remove()
        thetas = thetas_hist[frame]
        draw_compass_batch(ax, xs, ys, thetas,
                           length=needle_len, width=needle_width)
        draw_ext_field_on_lattice(ax, xs, ys, B_ext, phi_ext_deg, needle_len)
        ax.set_title(f"Relaxação dipolar — passo {frame * 20}",
                     color='#ECF0F1', fontsize=11, fontfamily='monospace')
        return []

    ani = FuncAnimation(fig, update, frames=len(thetas_hist),
                        interval=interval, blit=False)
    if save_gif:
        ani.save(save_gif, writer='pillow', fps=12, dpi=90)
        _print(f"GIF salvo em: {save_gif}")
    return fig, ani


# ══════════════════════════════════════════════════════════════════════════════
# 9. GLOBAL ORDER-PARAMETER PLOT
# ══════════════════════════════════════════════════════════════════════════════

def plot_order_parameter(thetas_hist, outpath, dt=None):
    """Plots the time evolution of the global magnetic order parameter S(t).

    S(t) = |<e^{i theta}>| in [0, 1]

S close to 1 means aligned needles; S close to 0 means random orientations.

Parameters
----------
thetas_hist : list of N×M arrays, the `relax` history sampled every 20 steps
outpath     : output-file path (PNG)
dt          : time step [s]; if provided, the x axis is in real seconds, otherwise it is the snapshot index"""
    order_params = [
        np.abs(np.mean(np.exp(1j * th)))
        for th in thetas_hist
    ]

    if dt is not None:
        # each snapshot corresponds to 20 steps
        time_ax = np.arange(len(order_params)) * 20 * dt
        xlabel  = "Tempo  t  [s]"
    else:
        time_ax = np.arange(len(order_params))
        xlabel  = "Snapshot (a cada 20 passos)"

    fig, ax = plt.subplots(figsize=(8, 3.5), facecolor='#1A1A2E')
    ax.set_facecolor('#0F3460')
    ax.plot(time_ax, order_params, color='#E94560', lw=2)
    ax.set_xlabel(xlabel, color='#BDC3C7')
    ax.set_ylabel(r"Parâmetro de ordem $S = |\langle e^{i\theta}\rangle|$",
                  color='#BDC3C7')
    ax.set_title("Evolução do parâmetro de ordem magnético global",
                 color='#ECF0F1', fontfamily='monospace')
    ax.tick_params(colors='#7F8C8D')
    ax.set_ylim(0, 1.05)
    ax.grid(True, color='#2C3E50', alpha=0.5)
    for spine in ax.spines.values():
        spine.set_edgecolor('#2C3E50')
    plt.tight_layout()
    plt.savefig(outpath, dpi=130, bbox_inches='tight', facecolor='#1A1A2E')
    _print(f"Parâmetro de ordem salvo em: {outpath}")


# ══════════════════════════════════════════════════════════════════════════════
# 10. INTERFACE OF LINE OF COMANDO (main)
# ══════════════════════════════════════════════════════════════════════════════

def main():
    """Command-line entry point of the simulation.

Execution flow
--------------
1. Parses arguments: geometry, spacing, external field, and related options
2. Converts (B_ext, phi_ext) into (Bext_x, Bext_y) Cartesian coordinates
3. Generates the grid of positions and initial angles with make_grid
4. Computes the cutoff and needle size according to geometry
5. Saves the initial-state figure
6. Runs relaxation with relax
7. Saves the equilibrium-state figure
8. Saves the initial-vs-equilibrium comparison figure
9. Plots the order parameter during relaxation
10. Optionally generates the relaxation animation"""
    parser = argparse.ArgumentParser(
        description="Simulacao de rede de bussolas — campo dipolar 2D",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # ── group 1: geometria of the lattice ────────────────────────────────────────
    grp = parser.add_argument_group("Geometria da rede")
    grp.add_argument('--N', type=int, default=8,
                     help='Numero de LINHAS de agulhas')
    grp.add_argument('--M', type=int, default=8,
                     help='Numero de COLUNAS de agulhas')
    grp.add_argument('--R', type=float, default=0.025,
                     help='Raio do circulo que envolve cada agulha [m]. '
                          'A distancia entre centros de vizinhos = 2R. '
                          'Padrao: 0.025 m = 2.5 cm')
    grp.add_argument('--needle_frac', type=float, default=0.80,
                     help='Comprimento da agulha como fracao do diametro 2R '
                          '(0.0 a 0.8). Padrao: 0.80  ->  agulha = 0.80 * 2R. '
                          'Valores fora do intervalo serao limitados (clamp).')
    grp.add_argument('--geometry',
                     choices=['square', 'triangular', 'honeycomb'],
                     default='square',
                     help='Tipo de rede')

    # ── group 2: physical and simulation parameters ────────────────────────
    grp2 = parser.add_argument_group("Fisica e simulacao")
    grp2.add_argument('--moment', type=float, default=None,
                      help='Momento magnetico de cada agulha [A·m²]. '
                           'Se omitido (padrao), calculado automaticamente '
                           'a partir do volume real da agulha (R, '
                           'needle_frac, espessura), assumindo que o aco '
                           'esta saturado na direcao do eixo longo (Norte) '
                           'da agulha: m = Ms * volume, onde Ms e a '
                           'magnetizacao de saturacao do aco (--steel_Bsat). '
                           'Se fornecido explicitamente, este valor '
                           'sobrescreve o calculo automatico. '
                           'Tipico: 0.01 (pequena) a 1.0 (grande).')
    grp2.add_argument('--inertia', type=float, default=None,
                      help='Momento de inercia de cada agulha [kg·m²]. '
                           'Se omitido (padrao), calculado automaticamente '
                           'a partir da geometria real da agulha (R, '
                           'needle_frac, espessura) como uma lamina de aco '
                           'magnetico: I = (1/12)*massa*(L²+largura²), '
                           'massa = densidade*L*largura*espessura. '
                           'Se fornecido explicitamente, este valor '
                           'sobrescreve o calculo automatico.')
    grp2.add_argument('--needle_thickness', type=float,
                      default=NEEDLE_THICKNESS_DEFAULT,
                      help='Espessura da lamina de aco da agulha [m], '
                           'usada no calculo automatico do momento de '
                           'inercia E do momento magnetico (quando '
                           '--inertia / --moment nao sao fornecidos). '
                           f'Padrao: {NEEDLE_THICKNESS_DEFAULT*1e3:.2f} mm')
    grp2.add_argument('--steel_density', type=float,
                      default=STEEL_DENSITY_DEFAULT,
                      help='Densidade do aco magnetico [kg/m³], usada '
                           'apenas no calculo automatico do momento de '
                           f'inercia. Padrao: {STEEL_DENSITY_DEFAULT:.0f} '
                           'kg/m³ (aco-carbono comum)')
    grp2.add_argument('--steel_Bsat', type=float,
                      default=None,
                      help='Densidade de fluxo de saturacao do aco '
                           'magnetico [T], usada apenas no calculo '
                           'automatico do momento magnetico (quando '
                           '--moment nao e fornecido). Relacao: '
                           'Ms = Bsat / mu0. Faixa tipica de acos '
                           'magneticos comuns: 1.6-2.2 T. '
                           f'Padrao: {STEEL_MS_SATURATION_DEFAULT*4*np.pi*1e-7:.2f} T')
    grp2.add_argument('--damping', type=float, default=DAMPING_DEFAULT,
                      help='Amortecimento viscoso do ar b [N·m·s/rad]. '
                           'Zero = sem amortecimento (oscila infinitamente). '
                           'Grande = sem oscilacoes visíveis. '
                           f'Padrao: {DAMPING_DEFAULT:.2e} N·m·s/rad '
                           '(sub-amortecido, Q >> 1)')
    grp2.add_argument('--t_sim', type=float, default=2.0,
                      help='Tempo fisico total da simulacao [s] '
                           '(= soma dos passos dt integrados). '
                           'Em hysteresis: cobre 1 ciclo completo (0->Hmax->-Hmax->Hmax). '
                           'Em sine: deve cobrir varios periodos (>= 3/field_freq). '
                           'Padrao: 2.0 s')
    grp2.add_argument('--t_sim_full', type=int, default=0, choices=[0, 1],
                      help='1 = roda ate o FIM de t_sim, desativando todas as '
                           'paradas antecipadas por criterio fisico (S=1.00, '
                           'rede em repouso, equilibrio por torque nos modos '
                           'pulse/step). Util para series temporais de '
                           'comprimento uniforme em varreduras. Interrupcoes '
                           'do usuario (Ctrl+I/Ctrl+C) continuam ativas. '
                           '0 (padrao) = paradas antecipadas normais.')
    grp2.add_argument('--field_mode',
                      choices=['static', 'hysteresis', 'sine', 'pulse',
                               'step_pos', 'step_neg'],
                      default='static',
                      help='Modo do campo externo: '
                           'static=constante (padrao), '
                           'hysteresis=ciclo completo 0->+Bmax->0->-Bmax->0->+Bmax (5 rampas), '
                           'sine=senoidal Hmax*sin(2pi*f*t), '
                           'pulse=espera field_delay, aplica por t_pulse (ou ate '
                           'S>=0.99 se t_pulse nao dado), zera e estabiliza, '
                           'step_pos=B=0 por field_delay depois campo LIGADO ate o fim, '
                           'step_neg=campo desde t=0, REMOVIDO em field_delay ate o fim. '
                           'Amplitude=--B_ext, direcao=--phi_ext.')
    grp2.add_argument('--hyst_spacing', choices=['linear', 'log'],
                      default='linear',
                      help='Espacamento da rampa de histerese: linear (padrao) '
                           'ou log (tipo log simetrico via sinh) -- dB/dt pequeno '
                           'perto de B=0 (amostragem FINA em campo na regiao da '
                           'transicao/coercividade) e grande perto de +-Bmax '
                           '(amostragem grossa na saturacao). '
                           'Apenas para --field_mode hysteresis.')
    grp2.add_argument('--hyst_log_k', type=float, default=5.0,
                      help='Concentracao do espacamento log da histerese: '
                           'k->0 recupera o linear; k=5 (padrao) concentra ~15x '
                           'mais pontos perto de B=0 do que na saturacao. '
                           'Apenas com --hyst_spacing log.')
    grp2.add_argument('--field_delay', type=float, default=0.0,
                      help='Delta_t de espera [s]: em pulse, tempo com B=0 antes '
                           'do pulso; em step_pos, tempo com B=0 antes de ligar; '
                           'em step_neg, tempo com campo LIGADO antes de remover. '
                           'Padrao: 0.')
    grp2.add_argument('--t_pulse', type=float, default=None,
                      help='Duracao do pulso [s] (apenas --field_mode pulse). '
                           'Se omitido, usa o criterio legado: campo ate S>=0.99 '
                           'ou plato de S.')
    grp2.add_argument('--torque_tol', type=float, default=1e-3,
                      help='Tolerancia relativa do criterio de equilibrio por '
                           'torque nos modos pulse/step_pos/step_neg: encerra '
                           'quando omega_max pequeno E media |tau| < tol * tau_ref '
                           '(tau_ref = I*omega0^2) por 50 passos. 0 desliga o '
                           'criterio de torque (mantem so omega). Padrao: 1e-3.')
    grp2.add_argument('--field_freq', type=float, default=1.0,
                      help='Frequencia do campo senoidal [Hz]. '
                           'Apenas para --field_mode sine. Padrao: 1.0 Hz')
    grp2.add_argument('--dt_factor', type=float, default=0.05,
                      help='Fracao do periodo natural T0 usada como passo dt '
                           '(0.02-0.10; menor = mais preciso, mais lento)')
    grp2.add_argument('--noise', type=float, default=1.5,
                      help='Amplitude do ruido inicial nos angulos [rad]; '
                           '0=todas para +x, 3.14=aleatorio total')
    grp2.add_argument('--seed', type=int, default=42,
                      help='Semente do gerador aleatorio (reprodutibilidade)')
    grp2.add_argument('--pbc', type=int, choices=[0, 1], default=0,
                      help='Condicoes periodicas de contorno (PBC): '
                           '0=desligado (padrao, rede finita com bordas), '
                           '1=ligado (rede tratada como estrutura periodica '
                           'infinita em x e y; agulhas nas bordas interagem '
                           'com replicas do lado oposto via convencao de '
                           'imagem minima)')
    grp2.add_argument('--pbc_images', type=int, default=1,
                      help='Numero de replicas periodicas somadas de cada '
                           'lado, em cada direcao, quando --pbc 1. Padrao: 1 '
                           '(soma sobre uma grade de (2*1+1)^2=9 celulas). '
                           'Usado apenas se --pbc 1.')
    grp2.add_argument('--gpu', type=int, choices=[0, 1], default=0,
                      help='Uso de GPU (CuPy) para o calculo do campo '
                           'dipolar: 0=desativado (padrao, forca CPU/NumPy '
                           'mesmo se uma GPU CUDA estiver disponivel), '
                           '1=ativado (usa GPU se CuPy estiver instalado e '
                           'funcional; cai para CPU automaticamente caso '
                           'contrario, com aviso no terminal)')
    grp2.add_argument('--progress_bar', type=int, choices=[0, 1], default=1,
                      help='Exibe a barra de progresso em tempo real durante '
                           'a integracao (ex: "Integrando [CPU] [###...] '
                           '50.0%% passo N/M"): 1=exibe (padrao), '
                           '0=nao exibe. Util para desabilitar em logs '
                           'redirecionados para arquivo, onde a atualizacao '
                           'por \\r nao se comporta como em um terminal.')
    grp2.add_argument('--halo_mode', type=str, choices=['order', 'domains'],
                      default='order',
                      help='Modo de coloracao dos halos em torno de cada '
                           'agulha: "order" (padrao) - verde/vermelho '
                           'conforme o parametro de ordem LOCAL (alinhamento '
                           'medio com vizinhos diretos); "domains" - cada '
                           'regiao CONEXA de agulhas com orientacao similar '
                           '(dentro de --domain_tol) recebe uma cor '
                           'categorica distinta, identificando dominios '
                           'magneticos separados por paredes de dominio.')
    grp2.add_argument('--domain_tol', type=float, default=15.0,
                      help='Tolerancia angular [graus] entre agulhas '
                           'vizinhas para serem consideradas parte do '
                           'mesmo dominio magnetico. Usado apenas quando '
                           '--halo_mode domains. Padrao: 15.0 graus')

    # ── group 3: external field ────────────────────────────────────────────
    grp3 = parser.add_argument_group(
        "Campo externo uniforme (SI)",
        "Intensidade em Tesla. Exemplos: campo terrestre = 50e-6 T; "
        "ima de geladeira a 5 cm ≈ 1e-3 T.")
    grp3.add_argument('--B_ext', type=float, default=0.0,
                      help='Intensidade do campo externo [T]. '
                           'Ex: 50e-6 (campo terrestre), 1e-3 (1 mT)')
    grp3.add_argument('--phi_ext', type=float, default=0.0,
                      help='Direcao do campo externo [graus]: '
                           '0=direita (+x), 90=cima (+y), anti-horario')
    grp3.add_argument('--ext_Bx', type=float, default=None,
                      help='Componente Bx do campo externo [T] '
                           '(sobrescreve --B_ext/--phi_ext se informado)')
    grp3.add_argument('--ext_By', type=float, default=None,
                      help='Componente By do campo externo [T] '
                           '(sobrescreve --B_ext/--phi_ext se informado)')

    # ── group 4: output ────────────────────────────────────────────────────
    grp4 = parser.add_argument_group("Saida")
    grp4.add_argument('--video', type=str, default=None,
                      metavar='NOME',
                      help='Gera video MP4 com este nome base (requer ffmpeg). '
                           'Extensao .mp4 e opcional e adicionada automaticamente. '
                           'Sempre recebe sufixo numerico: nome0000.mp4, '
                           'nome0001.mp4, etc. Ex: --video simples '
                           '-> simples0000.mp4')
    grp4.add_argument('--frame_every', type=int, default=5,
                      help='Salva um frame a cada N passos (menor = mais suave, '
                           'mais lento). Padrao: 5')
    grp4.add_argument('--make_images', type=int, default=1, choices=[0, 1],
                      help='1 (padrao) = gera as imagens PNG normalmente '
                           '(frames, compass_initial.png, compass_equilibrium.png, '
                           'compass_comparison.png, compass_order_param.png, e '
                           'hysteresis_loop.png/sine_field.png quando aplicavel). '
                           '0 = desliga TODA geracao de imagem/PNG, mantendo '
                           'apenas os CSVs de dados -- util para campanhas de '
                           'varredura de parametro com muitas execucoes, onde '
                           'o custo de renderizacao matplotlib domina o tempo '
                           'total e as imagens nao serao usadas. Nao afeta '
                           '--video: se --video for passado junto com '
                           '--make_images 0, nenhum frame existira para montar '
                           'o video e um aviso sera emitido.')
    grp4.add_argument('--fps', type=int, default=24,
                      help='Quadros por segundo do video MP4. Padrao: 24')
    grp4.add_argument('--dpi', type=int, default=120,
                      help='Resolucao dos frames em pontos por polegada. '
                           'Padrao: 120. Para redes grandes (>15x15) use '
                           '150-200 para melhor qualidade de imagem.')
    grp4.add_argument('--keep_frames', action='store_true',
                      help='Mantém a pasta de PNGs intermediários após gerar o MP4')
    grp4.add_argument('--csv_order', choices=['t', 'B'], default='t',
                      help='Ordem das colunas no CSV exportado: '
                           't (padrao) = tempo na 1a coluna: t,B,M_proj,S. '
                           'B = tempo na ultima coluna: B,M_proj,S,t.')
    args = parser.parse_args()

    # ── validates and limits needle_frac ao intervalo [0, 0.8] ──────────────────
    # above 0.8 the needles begin to touch or overlap visually
    if args.needle_frac < 0.0 or args.needle_frac > 0.8:
        clamped = max(0.0, min(args.needle_frac, 0.8))
        _print(f"  Aviso: --needle_frac {args.needle_frac} fora do intervalo "
                f"[0, 0.8]; ajustado para {clamped}")
        args.needle_frac = clamped

    np.random.seed(args.seed)

    # ── converts external field for coordeanythings cartesianas ───────────────
    # Prioridade: --ext_Bx/--ext_By (cartesiano) sobrescreve --B_ext/--phi_ext
    if args.ext_Bx is not None or args.ext_By is not None:
        Bext_x = args.ext_Bx if args.ext_Bx is not None else 0.0
        Bext_y = args.ext_By if args.ext_By is not None else 0.0
        B_ext_mag = np.sqrt(Bext_x**2 + Bext_y**2)
        phi_ext_deg = np.degrees(np.arctan2(Bext_y, Bext_x))
    else:
        # converts (intensidade, angle) → (Bx, By)
        phi_rad = np.deg2rad(args.phi_ext)
        Bext_x  = args.B_ext * np.cos(phi_rad)
        Bext_y  = args.B_ext * np.sin(phi_rad)
        B_ext_mag   = args.B_ext
        phi_ext_deg = args.phi_ext

    # ── sizes derived from R (computed here for use in the summary and in the
    #    automatic moment-of-inertia calculation, before printing) ──────
    # R  = radius of the circle visual of each needle
    # 2R = distance between neighbors (nn_dist) in any geometria
    # needle_len = fraction of 2R controlled by --needle_frac
    # needle_width = 22% of needle_len (visual compass proportion)
    R          = args.R
    needle_len = args.needle_frac * 2.0 * R   # length of the needle [m]
    needle_width = needle_len * 0.22           # width of the diamond [m]

    # ── moment of inertia: computed from geometry, or explicit value ────
    # Se the user not passou --inertia, computes automatically from
    # the real needle dimensions (length, width, thickness) and
    # the steel density, modeling the needle as a thin rectangular sheet.
    # If --inertia was provided explicitly, that value takes precedence.
    if args.inertia is None:
        inertia = compute_inertia_from_geometry(
            needle_len, needle_width, args.needle_thickness,
            density=args.steel_density)
        _inertia_auto = True
    else:
        inertia = args.inertia
        _inertia_auto = False

    # ── magnetic moment: computed from geometry (saturated steel), or
    #    explicit value ──────────────────────────────────────────────────
    # Se the user not passou --moment, computes automatically assumindo
    # that the needle is SATURATED along its long axis direction (North):
    # m = Ms * volume, with volume = L * width * thickness (same
    # geometry already used for the inertia). If --moment was provided
    # explicitly, that value takes precedence.
    if args.steel_Bsat is not None:
        _Ms_used = args.steel_Bsat / (4.0 * np.pi * 1e-7)   # Bsat = mu0 * Ms
    else:
        _Ms_used = STEEL_MS_SATURATION_DEFAULT
    if args.moment is None:
        moment = compute_moment_from_geometry(
            needle_len, needle_width, args.needle_thickness, Ms=_Ms_used)
        _moment_auto = True
    else:
        moment = args.moment
        _moment_auto = False

    # ── summary of the parameters ─────────────────────────────────────────────
    # formats field for display in a more readable unit
    def fmt_field(B):
        if B == 0:      return "0 T"
        if B >= 0.1:    return f"{B:.4f} T"
        if B >= 1e-4:   return f"{B*1e3:.4f} mT"
        return              f"{B*1e6:.2f} µT"

    _print(f"\n{'═'*62}")
    if args.gpu and _GPU_AVAILABLE:
        try:
            _gpu_name = cp.cuda.runtime.getDeviceProperties(0)['name'].decode()
        except Exception:
            _gpu_name = "GPU CUDA"
        _print(f"  Backend      : GPU ({_gpu_name}) via CuPy  [--gpu 1]")
    elif args.gpu and not _GPU_AVAILABLE:
        _print(f"  Backend      : CPU (NumPy)  [--gpu 1 solicitado, mas GPU indisponivel]")
        if _GPU_ERROR_MSG:
            if 'CUDA headers' in _GPU_ERROR_MSG or 'CUDA_PATH' in _GPU_ERROR_MSG:
                _print(f"  GPU indisponivel: faltam os headers do CUDA Toolkit.")
                _print(f"  Solucao: pip install cupy-cuda12x[ctk]")
                _print(f"  (ou ajuste cuda12x para a versao da sua CUDA Toolkit)")
            else:
                _print(f"  GPU indisponivel: {_GPU_ERROR_MSG[:70]}")
        else:
            _print(f"  GPU indisponivel: cupy nao instalado")
    else:
        _gpu_hint = " (GPU detectada, use --gpu 1 para ativar)" if _GPU_AVAILABLE else ""
        _print(f"  Backend      : CPU (NumPy)  [--gpu 0]{_gpu_hint}")
    _print(f"  Rede         : {args.geometry}  {args.N}x{args.M} agulhas")
    pbc_str = (f"ligado  (n_images={args.pbc_images}, soma sobre "
               f"{(2*args.pbc_images+1)**2} celulas)") if args.pbc else "desligado"
    _print(f"  PBC          : {pbc_str}")
    _print(f"  Raio R       : {args.R*100:.2f} cm  (2R = {2*args.R*100:.2f} cm)")
    _print(f"  Agulha       : {args.needle_frac*100:.0f}% de 2R = {args.needle_frac*2*args.R*100:.2f} cm")
    if _moment_auto:
        _Bsat_used = _Ms_used * 4.0 * np.pi * 1e-7
        _print(f"  Momento mag. : {moment:.4g} A.m2  por agulha  "
               f"[calculado: aco saturado a {_Bsat_used:.2f}T, "
               f"volume {needle_len*needle_width*args.needle_thickness*1e9:.3f}mm3]")
    else:
        _print(f"  Momento mag. : {moment:.4g} A.m2  por agulha  [valor manual]")
    if _inertia_auto:
        _mass_g = args.steel_density * needle_len * needle_width * args.needle_thickness * 1e3
        _print(f"  Inercia      : {inertia:.3e} kg.m2  por agulha  "
               f"[calculado: lamina de aco {args.needle_thickness*1e3:.2f}mm, "
               f"~{_mass_g:.3f}g]")
    else:
        _print(f"  Inercia      : {inertia:.3e} kg.m2  por agulha  [valor manual]")
    _print(f"  Amortecimento: {args.damping:.3e} N.m.s/rad  (ar)")
    _print(f"  Campo externo: {fmt_field(B_ext_mag)}  phi={phi_ext_deg:.1f} graus")
    _print(f"  Componentes  : Bx={fmt_field(abs(Bext_x))}  By={fmt_field(abs(Bext_y))}")
    field_mode_str = args.field_mode
    if args.field_mode == 'sine':
        field_mode_str = f"sine  f={args.field_freq:.3f} Hz"
    elif args.field_mode == 'hysteresis':
        field_mode_str = ("hysteresis  (0->+Bmax->0->-Bmax->0->+Bmax, "
                          f"espacamento {args.hyst_spacing}"
                          + (f" k={args.hyst_log_k:g}" if args.hyst_spacing == 'log' else "")
                          + ")")
    elif args.field_mode == 'pulse':
        if args.t_pulse is not None:
            field_mode_str = (f"pulse  (espera {args.field_delay:g}s -> campo por "
                              f"{args.t_pulse:g}s -> zera -> estabiliza)")
        else:
            field_mode_str = (f"pulse  (espera {args.field_delay:g}s -> campo ate "
                              f"S>=0.99/plato -> zera -> estabiliza)")
    elif args.field_mode == 'step_pos':
        field_mode_str = (f"step_pos  (B=0 por {args.field_delay:g}s -> campo "
                          f"LIGADO ate o fim)")
    elif args.field_mode == 'step_neg':
        field_mode_str = (f"step_neg  (campo desde t=0 -> REMOVIDO em "
                          f"{args.field_delay:g}s)")
    _print(f"  Modo campo   : {field_mode_str}")
    _print(f"  t_sim        : {args.t_sim:.3f} s  dt_factor={args.dt_factor}")
    _print(f"  Ruido/seed   : {args.noise:.2f} rad  seed={args.seed}")
    _print(f"{'═'*62}\n")

    # ── generates grid of positions and angles iniciais ─────────────────────────
    # make_grid uses R as the only parameter; also returns nn_dist = 2R
    # and Lx, Ly = lattice period for periodic boundary conditions
    xs, ys, thetas_init, nn_dist, Lx_period, Ly_period = make_grid(
        args.N, args.M,
        geometry=args.geometry,
        noise=args.noise,
        R=args.R,
    )

    # ── remaining sizes derived from R (R/needle_len/needle_width and the
    #    inertia were already computed before the terminal summary) ──────────
    # r_halo = R (exactly the circle radius) — circles touch tangentially
    # cutoff = 2.6 * 2R — cobre 1ª and 2ª layers of neighbors
    r_halo     = R * 0.98                      # slightly smaller than R for a visible gap
    cutoff     = nn_dist * 2.6                 # radius of cutoff of the interaction dipolar [m]

    # ── PBC: limits cutoff the min(Lx,Ly)/2 for evitar count dupla ───────
    # Above this limit, the minimum-image convention may count the same
    # replica more of a time, gerando physically incorrect results.
    if args.pbc:
        max_cutoff_pbc = min(Lx_period, Ly_period) / 2.0
        if cutoff > max_cutoff_pbc:
            _print(f"  PBC: cutoff reduzido de {cutoff*100:.2f}cm para "
                    f"{max_cutoff_pbc*100:.2f}cm (limite min(Lx,Ly)/2)")
            cutoff = max_cutoff_pbc

    # reference dipolar field between nearest neighbors (for display)
    B_ref = MU0_OVER_4PI * 2.0 * moment / nn_dist**3
    _print(f"  B_dipolar ref: {fmt_field(B_ref)}  (entre vizinhos)")
    if B_ext_mag > 0:
        ratio = B_ext_mag / B_ref
        dom = "DOMINANTE" if ratio > 1 else "fraco"
        _print(f"  B_ext/B_ref  : {ratio:.3f}  (campo externo {dom})")
    _print()

    ext_kwargs = dict(B_ext=B_ext_mag, phi_ext_deg=phi_ext_deg,
                      r_halo=r_halo, halo_mode=args.halo_mode,
                      domain_tol_deg=args.domain_tol)

    # ── figure size proportional to the lattice extent ─────────────────
    # Calculated once and reused in all frames and static figures.
    # The logic mirrors plot_state but uses xs/ys already known.
    _x_span = (xs.max() - xs.min()) + 4 * needle_len
    _y_span = (ys.max() - ys.min()) + 6 * needle_len
    _aspect = _x_span / _y_span if _y_span > 0 else 1.0
    _base   = 8.0
    if _aspect >= 1.0:
        _fw = min(_base * _aspect, 20.0)
        _fh = _fw / _aspect
    else:
        _fh = min(_base / _aspect, 20.0)
        _fw = _fh * _aspect
    figsize_main = (max(_fw, 6.0), max(_fh, 6.0))

    ext_kwargs['figsize_inches'] = figsize_main

    # ── figure of the estado initial ──────────────────────────────────────────
    if args.make_images:
        fig0, _ = plot_state(thetas_init, xs, ys,
                             title="Estado inicial (aleatório)",
                             needle_len=needle_len, needle_width=needle_width,
                             **ext_kwargs)
        plt.tight_layout()
        plt.savefig("compass_initial.png",
                    dpi=args.dpi, bbox_inches='tight', facecolor='#1A1A2E')
        plt.close(fig0)
        _print("Estado inicial salvo.")

    # ── inertial dynamic integration ──────────────────────────────────────
    frame_dir   = None
    final_video = None
    if args.video and not args.make_images:
        _print(f"  Aviso: --video '{args.video}' foi pedido junto com "
               f"--make_images 0. Nenhum frame sera gerado, logo nenhum "
               f"video sera montado. Ignorando --video.")
        args.video = None

    if args.video:
        import os
        # if the user did not provide an extension, automatically adds .mp4
        # ex: "--video simples" -> "simples.mp4" -> "simples0000.mp4"
        if not os.path.splitext(args.video)[1]:
            args.video = args.video + ".mp4"
        # resolves the final video name BEFORE creating the frames folder,
        # ensuring that the folder and file use the same base name
        final_video = next_available_path(args.video)
        _print(f"  Video sera salvo como '{final_video}'")
        base      = os.path.splitext(final_video)[0]
        frame_dir = base + "_frames"
        _print(f"Integrando e gravando frames em '{frame_dir}/'...")
    else:
        _print("Integrando dinâmica inercial...")

    thetas_eq, omegas_eq, thetas_hist, n_frames, sim_dt, stop_reason, field_log = relax(
        thetas_init.copy(), xs, ys,
        t_sim=args.t_sim,
        dt_factor=args.dt_factor,
        inertia=inertia,
        damping=args.damping,
        cutoff=cutoff,
        ext_field=(Bext_x, Bext_y),
        moment=moment,
        field_mode=args.field_mode,
        field_freq=args.field_freq,
        frame_dir=frame_dir,
        frame_every=args.frame_every,
        needle_len=needle_len,
        needle_width=needle_width,
        r_halo=r_halo,
        frame_dpi=args.dpi,
        figsize_inches=figsize_main,
        pbc=bool(args.pbc),
        Lx=Lx_period,
        Ly=Ly_period,
        n_images=args.pbc_images,
        B_ext=B_ext_mag,
        phi_ext_deg=phi_ext_deg,
        use_gpu=bool(args.gpu),
        show_progress=bool(args.progress_bar),
        halo_mode=args.halo_mode,
        domain_tol_deg=args.domain_tol,
        make_images=bool(args.make_images),
        hyst_spacing=args.hyst_spacing,
        hyst_log_k=args.hyst_log_k,
        field_delay=args.field_delay,
        t_pulse=args.t_pulse,
        torque_tol=args.torque_tol,
        t_sim_full=bool(args.t_sim_full),
    )
    frames_str = f"  ({n_frames} frames salvos)" if n_frames else ""
    _print(f"Integracao concluida - {stop_reason}{frames_str}")

    # ── exporta CSV universal with field and magnetization ────────────────────
    # Name = same name as the video (if any), otherwise "compass_field_log.csv"
    import csv as _csv, os as _os
    if args.video and final_video:
        csv_path = _os.path.splitext(final_video)[0] + ".csv"
    else:
        csv_path = "compass_field_log.csv"

    # ── order of the columns conforme --csv_order ────────────────────────────
    # field_log stores tuples (t, B, M_proj, S)
    # t: time in the first column  -> t, B, M_proj, S
    # B: time in the last column -> B, M_proj, S, t
    if args.csv_order == 'B':
        header = ['B_aplicado_T', 'M_proj', 'S', 't_s']
        rows   = [(B, Mp, S, t) for (t, B, Mp, S) in field_log]
    else:
        header = ['t_s', 'B_aplicado_T', 'M_proj', 'S']
        rows   = field_log

    with open(csv_path, 'w', newline='') as f:
        w = _csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    _print(f"  CSV salvo: {csv_path}  ({len(field_log)} pontos)")

    # ── figure of the estado of equilibrium ────────────────────────────────────
    if args.make_images:
        title_eq = (f"Equilíbrio — {args.geometry} {args.N}×{args.M}"
                    + (f"  |  B={fmt_field(B_ext_mag)} @ {phi_ext_deg:.0f}°"
                       if B_ext_mag > 0 else ""))
        fig1, _ = plot_state(thetas_eq, xs, ys,
                             title=title_eq,
                             needle_len=needle_len, needle_width=needle_width,
                             **ext_kwargs)
        plt.tight_layout()
        plt.savefig("compass_equilibrium.png",
                    dpi=args.dpi, bbox_inches='tight', facecolor='#1A1A2E')
        plt.close(fig1)
        _print("Estado de equilíbrio salvo.")

        # ── figure comparativa side the side ────────────────────────────────
        fig2, axes = plt.subplots(1, 2, figsize=(15, 7), facecolor='#1A1A2E')
        margin = needle_len * 1.6
        for ax_ in axes:
            ax_.set_facecolor('#16213E')
            ax_.set_aspect('equal')
            ax_.set_xlim(xs.min() - margin, xs.max() + margin)
            ax_.set_ylim(ys.min() - margin, ys.max() + margin * 2.5)
            ax_.tick_params(left=False, bottom=False,
                            labelleft=False, labelbottom=False)
            for sp in ax_.spines.values():
                sp.set_edgecolor('#2C3E50')

        axes[0].set_title("Estado inicial", color='#ECF0F1',
                          fontsize=12, fontfamily='monospace')
        axes[1].set_title("Equilíbrio dipolar", color='#ECF0F1',
                          fontsize=12, fontfamily='monospace')

        N_g, M_g = thetas_init.shape
        # ── halos (V43): same mode/tolerance as the other figures ──────
        # The comparison previously never drew halos in any mode; now it uses the
        # same plot_state routine — in domains mode, the title of each
        # panel ganha the count of domains daquele estado.
        _nd0, _ns0 = draw_halos_batch(axes[0], xs, ys, thetas_init,
                                      needle_len, r_halo=r_halo,
                                      halo_mode=args.halo_mode,
                                      domain_tol_deg=args.domain_tol)
        _nd1, _ns1 = draw_halos_batch(axes[1], xs, ys, thetas_eq,
                                      needle_len, r_halo=r_halo,
                                      halo_mode=args.halo_mode,
                                      domain_tol_deg=args.domain_tol)
        if args.halo_mode == 'domains':
            axes[0].set_title(
                f"Estado inicial — {_ns0} domínios principais "
                f"(de {_nd0})", color='#ECF0F1', fontsize=12,
                fontfamily='monospace')
            axes[1].set_title(
                f"Equilíbrio dipolar — {_ns1} domínios principais "
                f"(de {_nd1})", color='#ECF0F1', fontsize=12,
                fontfamily='monospace')

        draw_compass_batch(axes[0], xs, ys, thetas_init,
                           length=needle_len, width=needle_width)
        draw_compass_batch(axes[1], xs, ys, thetas_eq,
                           length=needle_len, width=needle_width)

        # external-field arrows in both panels
        draw_ext_field_on_lattice(axes[0], xs, ys, B_ext_mag, phi_ext_deg, needle_len)
        draw_ext_field_on_lattice(axes[1], xs, ys, B_ext_mag, phi_ext_deg, needle_len)

        bfield_str = (f"  |  B_ext={fmt_field(B_ext_mag)} @ {phi_ext_deg:.0f}°"
                      if B_ext_mag > 0 else "")
        fig2.suptitle(
            f"Rede {args.geometry} {args.N}×{args.M}"
            f"  |  R={args.R*100:.1f} cm  2R={2*args.R*100:.1f} cm"
            f"{bfield_str}  |  interação dipolar 2D",
            color='#BDC3C7', fontsize=11, fontfamily='monospace')
        plt.tight_layout()
        plt.savefig("compass_comparison.png",
                    dpi=args.dpi, bbox_inches='tight', facecolor='#1A1A2E')
        plt.close(fig2)
        _print("Comparação salva.")

        # ── order parameter global ────────────────────────────────────
        # hist contains tuples (thetas, omegas); we extract only the thetas
        plot_order_parameter([th for th, _ in thetas_hist],
                             "compass_order_param.png",
                             dt=sim_dt)

    # ── MP4 video generation ──────────────────────────────────────────────
    if final_video and frame_dir and n_frames > 0:
        import shutil
        ok = render_video(frame_dir, final_video, fps=args.fps, use_gpu=bool(args.gpu))
        if ok and not args.keep_frames:
            shutil.rmtree(frame_dir)
            _print(f"  Pasta de frames removida: {frame_dir}/")
        elif not ok:
            _print(f"  Os frames PNG foram mantidos em: {frame_dir}/")

    _print("\nConcluído.")


if __name__ == '__main__':
    main()
