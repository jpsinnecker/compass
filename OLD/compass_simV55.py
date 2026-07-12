"""
============================================================
compass.py — Simulation of a compass-needle lattice V55
============================================================

Models a 2D grid of classical magnetic dipoles, represented as compass
needles, that interact through the magnetic field each needle produces at
its neighbors. The dynamics are inertial, using Newton's second law for
rotation with no pivot friction and viscous air damping. The integrator is
Velocity-Verlet. All physical quantities use SI units.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
COMMAND-LINE PARAMETERS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

┌─ LATTICE GEOMETRY ─────────────────────────────────────────┐
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

┌─ PHYSICS AND SIMULATION ───────────────────────────────────┐
│                                                            │
│  --moment    float Magnetic moment of each needle [A·m²]   │
│                    Default: 0.1 (table compass, about 5 cm)│
│                    Ref: pocket ≈ 0.01 | nautical ≈ 1.0     │
│                                                            │
│  --inertia   float Moment of inertia [kg·m²]. If omitted,  │
│                    it is computed automatically from the   │
│                    geometry: steel sheet using R,          │
│                    needle_frac, --needle_thickness,        │
│                    and --steel_density.                    │
│                                                            │
│  --damping   float Viscous air damping [N·m·s/rad]         │
│                    Controls the quality factor Q:          │
│                    Q = omega_0·I/b (high Q = more          │
│                    oscillatory motion). Default: 5e-8      │
│                    (Q≈25, realistic compass). For a smooth │
│                    B_ext=0.1 T run, use 8e-6 (Q≈4).        │
│                                                            │
│  --t_sim     float Total physical simulation time [s].     │
│                    Sum of all integrated dt steps,         │
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
│  --ext_Bx    float x component of the field [T]            │
│  --ext_By    float y component of the field [T]            │
└────────────────────────────────────────────────────────────┘

┌─ OUTPUT ───────────────────────────────────────────────────┐
│  PNG files are always generated in the current directory:  │
│    compass_initial.png      initial state                  │
│    compass_equilibrium.png  final state                    │
│    compass_comparison.png   side-by-side comparison        │
│    compass_order_param.png  order parameter S(t)           │
│                                                            │
│  --video     str   Path of the MP4 video to generate.      │
│                    Requires ffmpeg. If the file already    │
│                    exists, saves as name0001.mp4,          │
│                    name0002.mp4, etc.                      │
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

┌─ CONTROLS DURING THE SIMULATION ───────────────────────────┐
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
  python compass.py

  # Honeycomb with Earth's field
  python compass.py --geometry honeycomb --N 10 --M 10 --B_ext 50e-6

  # Triangular lattice with a 0.1 T field at 45°, smooth motion
  python compass.py --geometry triangular --N 10 --M 10 \
      --B_ext 0.1 --phi_ext 45 --damping 8e-6 --t_sim 2.0

  # Video with larger needles and many oscillations
  python compass.py --R 0.03 --needle_frac 0.85 --damping 1e-9 \
      --t_sim 5.0 --frame_every 2 --fps 30 --video sim.mp4

  # Field through Cartesian components
  python compass.py --ext_Bx 0.05 --ext_By -0.05

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

_GPU_AVAILABLE = False
_GPU_ERROR_MSG = None
_xp = np
try:
    import cupy as cp
    cp.cuda.Device(0).compute_capability
    _test_arr = cp.array([0.0, 1.0])
    _ = cp.cos(_test_arr)
    cp.cuda.Stream.null.synchronize()
    _xp = cp
    _GPU_AVAILABLE = True
except Exception as _e:
    _xp = np
    _GPU_AVAILABLE = False
    _GPU_ERROR_MSG = str(_e).strip().splitlines()[-1] if str(_e).strip() else type(_e).__name__


def _to_cpu(arr):
    if _GPU_AVAILABLE and hasattr(arr, 'get'):
        return arr.get()
    return arr


def _to_backend(arr):
    if _GPU_AVAILABLE:
        return cp.asarray(arr)
    return np.asarray(arr)


sys.stdout.reconfigure(line_buffering=True)

def _print(*args, **kwargs):
    import sys
    msg = ' '.join(str(a) for a in args)
    sys.stdout.write(msg + '\r\n')
    sys.stdout.flush()


def _progress_ansi_ok():
    import sys
    if _PROGRESS_REGION['ansi'] is None:
        try:
            _PROGRESS_REGION['ansi'] = bool(sys.stdout.isatty())
        except Exception:
            _PROGRESS_REGION['ansi'] = False
    return _PROGRESS_REGION['ansi']


_PROGRESS_REGION = {'open': False, 'ansi': None}


def _print_progress_bar(frac, bar_width=40, step=0, n_steps=0, steps_per_s=0,
                        t_now=0.0, B_status=0.0, S_now=0.0, omega_max=0.0,
                        backend="CPU"):
    import sys
    frac = max(0.0, min(1.0, frac))
    n_filled = int(round(frac * bar_width))
    bar = '#' * n_filled + ' ' * (bar_width - n_filled)
    pct = frac * 100.0
    
    line1 = f"Integrating [{backend}] [{bar}]"
    line2 = f"{pct:5.1f}%  step {step}/{n_steps}  ({steps_per_s:.0f} steps/s)  t={t_now:.4f}s"
    line3 = f"B={B_status*1e3:+.4f}mT  S={S_now:.4f}  w_max={omega_max:.3f}rad/s"

    if _progress_ansi_ok():
        seq = ("\x1b[K" + line1 + "\n\r" +
               "\x1b[K" + line2 + "\n\r" +
               "\x1b[K" + line3 + "\x1b[2A\r")
        sys.stdout.write(seq)
        sys.stdout.flush()
        _PROGRESS_REGION['open'] = True
    else:
        sys.stdout.write("\r" + line1 + " | " + line2 + " | " + line3)
        sys.stdout.flush()


def _print_progress_bar_finish():
    import sys
    if _progress_ansi_ok():
        if _PROGRESS_REGION['open']:
            sys.stdout.write('\n\r\n\r\n\r')
            sys.stdout.flush()
            _PROGRESS_REGION['open'] = False
    else:
        sys.stdout.write('\r\n')
        sys.stdout.flush()

# ══════════════════════════════════════════════════════════════════════════════
# 1. PHYSICAL CONSTANTS (International System — SI)
# ══════════════════════════════════════════════════════════════════════════════

MU0_OVER_4PI = 1.0e-7   # T·m/A
MOMENT_DEFAULT = 0.1      # A·m²
INERTIA_DEFAULT = 1.0e-7  # kg·m²
STEEL_DENSITY_DEFAULT = 7850.0   # kg/m³
NEEDLE_THICKNESS_DEFAULT = 0.4e-3   # m
STEEL_MS_SATURATION_DEFAULT = 1.59e6   # A/m


def compute_inertia_from_geometry(needle_len, needle_width, thickness,
                                  density=STEEL_DENSITY_DEFAULT):
    mass = density * needle_len * needle_width * thickness
    I = (1.0 / 12.0) * mass * (needle_len**2 + needle_width**2)
    return I


def compute_moment_from_geometry(needle_len, needle_width, thickness,
                                 Ms=STEEL_MS_SATURATION_DEFAULT):
    volume = needle_len * needle_width * thickness
    m = Ms * volume
    return m

DAMPING_DEFAULT = 5.0e-8  # N·m·s/rad


# ══════════════════════════════════════════════════════════════════════════════
# 2. MAGNETIC DIPOLAR FIELD (SI)
# ══════════════════════════════════════════════════════════════════════════════

def dipole_field_2d(rx, ry, theta_src, moment):
    r2 = rx**2 + ry**2
    if r2 < 1e-24:
        return 0.0, 0.0

    r  = np.sqrt(r2)
    r5 = r2 * r2 * r

    mx = moment * np.cos(theta_src)
    my = moment * np.sin(theta_src)

    mdotr = mx * rx + my * ry

    Bx = MU0_OVER_4PI * (3.0 * mdotr * rx / r5  -  mx / (r2 * r))
    By = MU0_OVER_4PI * (3.0 * mdotr * ry / r5  -  my / (r2 * r))
    return Bx, By


# ══════════════════════════════════════════════════════════════════════════════
# 3. TOTAL FIELD ON THE NEEDLE (SI)
# ══════════════════════════════════════════════════════════════════════════════

def total_field_on(i, j, thetas, xs, ys, cutoff, moment,
                   pbc=False, Lx=None, Ly=None, n_images=1):
    Bx_tot, By_tot = 0.0, 0.0
    N, M = thetas.shape
    xi, yi = xs[i, j], ys[i, j]

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
                continue

            rx0 = xi - xs[ni, nj]
            ry0 = yi - ys[ni, nj]

            for dx_img in x_shifts:
                for dy_img in y_shifts:
                    if ni == i and nj == j and dx_img == 0.0 and dy_img == 0.0:
                        continue

                    rx = rx0 + dx_img
                    ry = ry0 + dy_img
                    dist = np.sqrt(rx*rx + ry*ry)

                    if dist > cutoff:
                        continue

                    bx, by = dipole_field_2d(rx, ry, thetas[ni, nj], moment)
                    Bx_tot += bx
                    By_tot += by

    return Bx_tot, By_tot


# ══════════════════════════════════════════════════════════════════════════════
# 4. INERTIAL DYNAMICS (Newton's second law for rotation — no pivot friction)
# ══════════════════════════════════════════════════════════════════════════════

def _plot_hysteresis(log, hyst_autoscale=False):
    log = np.array(log)
    t_arr  = log[:, 0]
    B_arr  = log[:, 1] * 1e3
    M_arr  = log[:, 2]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), facecolor='#1A1A2E')

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

    if not hyst_autoscale:
        B_max = np.max(np.abs(B_arr))
        if B_max < 1e-5:
            B_max = 1.0
        ax1.set_xlim(-1.05 * B_max, 1.05 * B_max)
        ax1.set_ylim(-1.15, 1.15)

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

    fig.suptitle("Magnetic Hysteresis Simulation — Compass Array",
                 color='#BDC3C7', fontsize=12, fontfamily='monospace')
    plt.tight_layout()
    plt.savefig("hysteresis_loop.png", dpi=130, bbox_inches='tight',
                facecolor='#1A1A2E')
    plt.close(fig)
    _print("  Hysteresis plot saved: hysteresis_loop.png")


def _plot_sine(log, freq):
    log   = np.array(log)
    t_arr = log[:, 0]
    B_arr = log[:, 1] * 1e3
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
    _print("  Sinusoidal plot saved: sine_field.png")


# ══════════════════════════════════════════════════════════════════════════════
# 3b. VECTORIZED CALCULATION OF THE DIPOLAR FIELD AND TORQUES (GPU/CPU)
# ══════════════════════════════════════════════════════════════════════════════

def compute_torques_vectorized(theta_flat, x_flat, y_flat, moment, cutoff,
                               bx_ext, by_ext,
                               pbc=False, Lx=None, Ly=None, n_images=1):
    if _GPU_AVAILABLE and hasattr(theta_flat, 'get'):
        xp = cp
    else:
        xp = np
    K  = theta_flat.shape[0]

    mx = moment * xp.cos(theta_flat)
    my = moment * xp.sin(theta_flat)

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
            rx = (x_flat[:, None] - x_flat[None, :]) - dx_img
            ry = (y_flat[:, None] - y_flat[None, :]) - dy_img

            r2 = rx*rx + ry*ry

            is_self = (r2 < 1e-24)
            valid   = (~is_self) & (r2 <= cutoff*cutoff)

            r2_safe = xp.where(valid, r2, 1.0)
            r       = xp.sqrt(r2_safe)
            r5      = r2_safe * r2_safe * r

            mdotr = mx[None, :] * rx + my[None, :] * ry

            bx_pair = MU0_OVER_4PI * (3.0 * mdotr * rx / r5 - mx[None, :] / (r2_safe * r))
            by_pair = MU0_OVER_4PI * (3.0 * mdotr * ry / r5 - my[None, :] / (r2_safe * r))

            bx_pair = xp.where(valid, bx_pair, 0.0)
            by_pair = xp.where(valid, by_pair, 0.0)

            Bx_tot += bx_pair.sum(axis=1)
            By_tot += by_pair.sum(axis=1)

    Bx_tot += bx_ext
    By_tot += by_ext

    tau_flat = mx * By_tot - my * Bx_tot
    return tau_flat


# ── V40: precomputed dipolar-interaction tensor ─────────────────────────
TENSOR_MEM_LIMIT_BYTES = 4.0e9


def precompute_dipolar_tensor(x_flat, y_flat, cutoff,
                              pbc=False, Lx=None, Ly=None, n_images=1):
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
          t_sim_full=False, hyst_autoscale=False):
    import os
    import time as _time_module

    N, M = thetas.shape
    Bext_x, Bext_y = ext_field
    n_frames = 0

    _active_gpu = use_gpu and _GPU_AVAILABLE
    _active_xp  = cp if _active_gpu else np

    def _local_to_backend(arr):
        if _active_gpu:
            return cp.asarray(arr)
        return np.asarray(arr)

    def _local_to_cpu(arr):
        if _active_gpu and hasattr(arr, 'get'):
            return arr.get()
        return arr

    x_flat_for_nn = xs.ravel()
    y_flat_for_nn = ys.ravel()
    dx_nn = x_flat_for_nn[:, None] - x_flat_for_nn[None, :]
    dy_nn = y_flat_for_nn[:, None] - y_flat_for_nn[None, :]
    d_nn  = np.sqrt(dx_nn**2 + dy_nn**2)
    np.fill_diagonal(d_nn, np.inf)
    r_nn  = float(np.min(d_nn))

    B_ref     = MU0_OVER_4PI * 2.0 * moment / r_nn**3
    B_ext_mag = np.sqrt(Bext_x**2 + Bext_y**2)
    B_eff     = max(B_ref, B_ext_mag)
    omega0    = np.sqrt(moment * B_eff / inertia)
    T0        = 2.0 * np.pi / omega0
    dt        = dt_factor * T0
    n_steps   = max(1, int(np.ceil(t_sim / dt)))
    Q         = omega0 * inertia / damping if damping > 0 else np.inf

    _print(f"  Dynamical inertial parameters:")
    _print(f"    r_nn    = {r_nn*100:.2f} cm")
    _print(f"    B_ref   = {B_ref*1e3:.4f} mT  (dipolar neighbours)")
    if B_ext_mag > 0:
        _print(f"    B_ext   = {B_ext_mag*1e3:.4f} mT  (External Field)")
    _print(f"    B_eff   = {B_eff*1e3:.4f} mT  (Dominant Field)")
    _print(f"    omega_0 = {omega0:.2f} rad/s   T0 = {T0:.5f} s")
    _print(f"    dt      = {dt:.6f} s  ({dt_factor:.0%} of T0)")
    _print(f"    t_sim   = {t_sim:.3f} s  -> {n_steps} steps")
    if Q > 2:
        q_desc = "sub-damping(oscilates)"
    elif Q > 0.5:
        q_desc = "critically damped"
    else:
        q_desc = "overdamped"
    _print(f"    Q       = {Q:.1f}  ({q_desc})")
    _print()

    theta_cur = thetas.copy()
    omega_cur = np.zeros((N, M))
    hist      = [(theta_cur.copy(), omega_cur.copy())]

    if not make_images and frame_dir is not None:
        _print(f"  Aviso: make_images=False -> ignorando frame_dir='{frame_dir}' "
               f"(nenhum PNG sera gerado)")
        frame_dir = None

    if frame_dir is not None:
        os.makedirs(frame_dir, exist_ok=True)

    def _fmt_B(B):
        if B == 0:    return ""
        if B >= 0.1:  return f"B={B:.3f} T"
        if B >= 1e-4: return f"B={B*1e3:.3f} mT"
        return            f"B={B*1e6:.1f} µT"

    def _draw_clock(ax, t_phys, needle_len, stop_label=None):
        xlim  = ax.get_xlim()
        ylim  = ax.get_ylim()
        xspan = xlim[1] - xlim[0]
        yspan = ylim[1] - ylim[0]

        px    = xlim[0] + 0.02 * xspan
        py_top = ylim[1] - 0.035 * yspan

        pad = needle_len * 0.25
        box_h = 0.140 * yspan + pad
        bar_w = 0.30 * xspan
        bar_h = 0.018 * yspan
        bar_y = py_top - 0.115 * yspan

        ax.add_patch(plt.Rectangle(
            (px - pad * 0.3, py_top - box_h),
            bar_w + pad, box_h,
            facecolor='#080818', edgecolor='#3A3A6A',
            linewidth=0.8, alpha=0.82, zorder=19,
            transform=ax.transData))

        frac = min(t_phys / t_sim, 1.0) if t_sim > 0 else 0.0

        ax.add_patch(plt.Rectangle(
            (px, bar_y), bar_w, bar_h,
            facecolor='#252540', edgecolor='none',
            zorder=20, transform=ax.transData))

        if frac > 0:
            r_col = min(2.0 * frac, 1.0)
            g_col = min(2.0 * (1.0 - frac), 1.0)
            ax.add_patch(plt.Rectangle(
                (px, bar_y), bar_w * frac, bar_h,
                facecolor=(r_col, g_col, 0.15), edgecolor='none',
                zorder=21, transform=ax.transData))

        ax.add_patch(plt.Rectangle(
            (px, bar_y), bar_w, bar_h,
            facecolor='none', edgecolor='#5A5A8A',
            linewidth=0.7, zorder=22, transform=ax.transData))

        if 0 < frac < 1.0:
            mx = px + bar_w * frac
            ax.plot([mx, mx], [bar_y, bar_y + bar_h],
                    color='white', lw=1.2, zorder=24,
                    transform=ax.transData)

        if t_phys >= 60.0:
            mins = int(t_phys // 60)
            secs = t_phys - mins * 60
            time_str = f"t = {mins:02d}:{secs:05.2f}"
        else:
            time_str = f"t = {t_phys:.4f} s"

        ax.text(px + pad * 0.2, py_top - 0.015 * yspan,
                time_str,
                color='#E8E8FF', fontsize=11, fontweight='bold',
                fontfamily='monospace', va='top', ha='left',
                zorder=25, transform=ax.transData)

        ax.text(px + pad * 0.2, py_top - 0.065 * yspan,
                f"/ {t_sim:.4f} s  (tempo físico)",
                color='#5555AA', fontsize=6,
                fontfamily='monospace', va='top', ha='left',
                zorder=25, transform=ax.transData)

        if stop_label:
            ax.text(px + pad * 0.2, py_top - 0.090 * yspan,
                    stop_label,
                    color='#FFD700', fontsize=7, fontweight='bold',
                    fontfamily='monospace', va='top', ha='left',
                    zorder=25, transform=ax.transData)

    field_log = []
    hysteresis_log = [] if field_mode == 'hysteresis' else None
    sine_log       = [] if field_mode == 'sine'       else None

    def _save_frame(step, th, om, stop_label=None,
                    B_ext_inst=None, phi_ext_inst=None):
        nonlocal n_frames
        t_phys   = step * dt
        S        = np.abs(np.mean(np.exp(1j * th)))
        om_max   = np.max(np.abs(om))

        b_now    = B_ext_inst   if B_ext_inst   is not None else B_ext
        phi_now  = phi_ext_inst if phi_ext_inst is not None else phi_ext_deg

        if b_now is not None and abs(b_now) > 1e-12:
            sign_str = "+" if b_now > 0 else "-"
            b_str = sign_str + _fmt_B(abs(b_now))
        else:
            b_str = "0 T"

        title = f"S = {S:.4f}   w_max = {om_max:.2f} rad/s   B = {b_str}"

        if field_mode == 'hysteresis':
            fig, (axL, axR) = plt.subplots(
                1, 2, figsize=(figsize_inches[0] * 2.05, figsize_inches[1]),
                facecolor='#1A1A2E',
                gridspec_kw={'width_ratios': [1.0, 1.0]})

            axL.set_facecolor('#16213E')
            axL.set_aspect('equal')
            _margin = needle_len * 1.6
            yspan_lattice = ys.max() - ys.min() + 2 * _margin
            yspan = yspan_lattice / 0.8
            axL.set_xlim(xs.min() - _margin, xs.max() + _margin)
            axL.set_ylim(ys.min() - _margin, ys.min() - _margin + yspan)
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
                phi_now, needle_len, B_ext_max=B_ext, B_signed=b_now)
            
            axL.set_title(title, color='#ECF0F1', fontsize=11,
                          fontfamily='monospace')
            _draw_clock(axL, t_phys, needle_len, stop_label=stop_label)

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
                             phi_ext_deg=phi_now,
                             B_ext_max=B_ext,
                             B_signed=b_now,
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
        axR.set_facecolor('#16213E')
        for sp in axR.spines.values():
            sp.set_edgecolor('#2C3E50')
        axR.tick_params(colors='#7F8C9A', labelsize=8)

        H_max_mT = (B_ext * 1e3) if B_ext > 0 else 1.0
        axR.set_xlabel("H  (mT)", color='#BDC3C7', fontsize=10,
                       fontfamily='monospace')
        axR.set_ylabel("M  (projeção normalizada)", color='#BDC3C7',
                       fontsize=10, fontfamily='monospace')
        axR.set_title("Curva de histerese  M × H", color='#ECF0F1',
                       fontsize=11, fontfamily='monospace')

        axR.axhline(0.0, color='#2C3E50', lw=1.0, zorder=1)
        axR.axvline(0.0, color='#2C3E50', lw=1.0, zorder=1)

        if hysteresis_log:
            arr = np.array(hysteresis_log)
            H_mT = arr[:, 1] * 1e3
            Mv   = arr[:, 2]
            if hyst_autoscale:
                H_min, H_max = np.min(H_mT), np.max(H_mT)
                M_min, M_max = np.min(Mv), np.max(Mv)
                H_span = H_max - H_min
                M_span = M_max - M_min
                
                H_margin = 0.05 * H_span if H_span > 1e-5 else 1.0
                M_margin = 0.05 * M_span if M_span > 1e-5 else 0.1
                
                axR.set_xlim(H_min - H_margin, H_max + H_margin)
                axR.set_ylim(M_min - M_margin, M_max + M_margin)
            else:
                axR.set_xlim(-1.05 * H_max_mT, 1.05 * H_max_mT)
                axR.set_ylim(-1.15, 1.15)
                
            axR.plot(H_mT, Mv, color='#F1C40F', lw=1.6, zorder=3)
            axR.plot([H_mT[-1]], [Mv[-1]], 'o', color='#E74C3C',
                     ms=7, zorder=4)
        else:
            axR.set_xlim(-1.05 * H_max_mT, 1.05 * H_max_mT)
            axR.set_ylim(-1.15, 1.15)

    if frame_dir is not None:
        _save_frame(0, theta_cur, omega_cur)

    import threading, sys, os as _os

    _stop_flag = threading.Event()
    _term_fd  = None
    _term_old = None
    try:
        import termios
        _term_fd  = sys.stdin.fileno()
        _term_old = termios.tcgetattr(_term_fd)
    except Exception:
        pass

    def _restore_terminal():
        if _term_fd is not None and _term_old is not None:
            try:
                import termios
                termios.tcsetattr(_term_fd, termios.TCSADRAIN, _term_old)
            except Exception:
                pass

    import atexit
    atexit.register(_restore_terminal)

    def _keyboard_listener():
        try:
            import tty, termios
            fd = sys.stdin.fileno()
            tty.setraw(fd)
            while not _stop_flag.is_set():
                ch = _os.read(fd, 1)
                if ch == b'\t':
                    _stop_flag.set()
                    break
                elif ch == b'\x03':
                    _restore_terminal()
                    sys.stdout.write("\r\n  Aborted (Ctrl+C) — immediate shutdown.\r\n")
                    sys.stdout.flush()
                    _os._exit(130)
        except Exception:
            pass
        finally:
            _restore_terminal()

    _kb_thread = threading.Thread(target=_keyboard_listener, daemon=True)
    _kb_thread.start()

    Bext_x0, Bext_y0 = ext_field
    B_max = np.sqrt(Bext_x0**2 + Bext_y0**2)
    phi_rad = np.arctan2(Bext_y0, Bext_x0)
    cos_phi = np.cos(phi_rad)
    sin_phi = np.sin(phi_rad)

    def field_at(t):
        if field_mode == 'static':
            return Bext_x0, Bext_y0

        elif field_mode == 'hysteresis':
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
            if t < field_delay:
                return 0.0, 0.0
            return Bext_x0, Bext_y0

        elif field_mode == 'step_neg':
            if t < field_delay:
                return Bext_x0, Bext_y0
            return 0.0, 0.0

        else:
            return Bext_x0, Bext_y0

    _pulse_relaxing = [False]

    _print("  [Ctrl+I to stop and e save video]")
    _print("  [Ctrl+C to ABORT withou saving]")

    x_flat_xp = _local_to_backend(xs.ravel())
    y_flat_xp = _local_to_backend(ys.ravel())

    _dipolar_tensor = None
    _K_total = int(x_flat_xp.shape[0])
    _tensor_bytes = 3.0 * _K_total * _K_total * 8.0
    if _tensor_bytes <= TENSOR_MEM_LIMIT_BYTES:
        _dipolar_tensor = precompute_dipolar_tensor(
            x_flat_xp, y_flat_xp, cutoff,
            pbc=pbc, Lx=Lx, Ly=Ly, n_images=n_images)
        _print(f"  Pre calculated Diploar Tensor: 3 matrixes {_K_total}x{_K_total} "
               f"({_tensor_bytes/1e6:.0f} MB) — step matriz-vetor (V40)")
    else:
        _print(f"  Dipolar tensor NOT precomputed ({_tensor_bytes/1e9:.1f} GB "
               f"exceeds the limit of {TENSOR_MEM_LIMIT_BYTES/1e9:.0f} GB) — "
               f"using step-by-step calculation (V39 method)")

    def _torques_xp(theta_flat_xp, bx_ext, by_ext):
        if _dipolar_tensor is not None:
            return compute_torques_from_tensor(
                theta_flat_xp, _dipolar_tensor, moment, bx_ext, by_ext)
        return compute_torques_vectorized(
            theta_flat_xp, x_flat_xp, y_flat_xp, moment, cutoff,
            bx_ext, by_ext, pbc=pbc, Lx=Lx, Ly=Ly, n_images=n_images)

    def _torques(th, bx_ext, by_ext):
        theta_flat_xp = _local_to_backend(th.ravel())
        tau_flat_xp = _torques_xp(theta_flat_xp, bx_ext, by_ext)
        tau_flat = _local_to_cpu(tau_flat_xp)
        return tau_flat.reshape(N, M)

    if _active_gpu:
        cp.cuda.Stream.null.synchronize()
    _perf_t_start  = _time_module.perf_counter()
    _perf_last_print_step = 0
    _perf_last_print_time = _perf_t_start
    _last_status_print_time = _perf_t_start
    _status_text = f"  [step 0/{n_steps}]  waiting for first status..."
    _last_bar_suffix = f"step 0/{n_steps}"

    def _finish_progress_bar_if_shown():
        if show_progress:
            _print_progress_bar_finish()

    bx_cur, by_cur    = field_at(0.0)
    theta_xp          = _local_to_backend(theta_cur.ravel())
    omega_xp           = _active_xp.zeros(N * M)
    tau_xp             = _torques_xp(theta_xp, bx_cur, by_cur)
    _converged_count  = 0
    _S1_count         = 0
    _allow_S1_stop     = (field_mode == 'static') and (not t_sim_full)
    _allow_rest_stop   = (field_mode == 'static') and (not t_sim_full)
    _stop_reason       = "tempo total atingido"
    _pulse_phase   = "delay" if (field_mode == 'pulse' and field_delay > 0) \
                     else "field_on"
    _S99_count        = 0
    _S_window         = []
    _tau_ref = inertia * omega0 * omega0
    _step_after_announced = False
    _S_WINDOW_SIZE    = max(40, int(2.0 / dt_factor))

    for step in range(1, n_steps + 1):

        t_now = step * dt

        if _stop_flag.is_set():
            _finish_progress_bar_if_shown()
            _print(f"\n  Stopped (Ctrl+I) at t={t_now:.4f}s  (step {step}/{n_steps})")
            _stop_reason = "interrupted by user (Ctrl+I)"
            if frame_dir is not None:
                theta_cur = _local_to_cpu(theta_xp).reshape(N, M)
                omega_cur = _local_to_cpu(omega_xp).reshape(N, M)
                _save_frame(step, theta_cur, omega_cur,
                            stop_label="■ interrupted (Ctrl+I)")
            break

        bx_new, by_new = field_at(t_now)

        accel_xp  = (tau_xp - damping * omega_xp) / inertia
        theta_new_xp = theta_xp + omega_xp * dt + 0.5 * accel_xp * dt**2
        theta_new_xp = (theta_new_xp + _active_xp.pi) % (2.0 * _active_xp.pi) - _active_xp.pi

        tau_new_xp = _torques_xp(theta_new_xp, bx_new, by_new)

        b_half    = damping * dt / (2.0 * inertia)
        omega_new_xp = (omega_xp * (1.0 - b_half)
                        + dt * (tau_xp + tau_new_xp) / (2.0 * inertia)) \
                       / (1.0 + b_half)

        theta_xp = theta_new_xp
        omega_xp = omega_new_xp
        tau_xp   = tau_new_xp
        bx_cur    = bx_new
        by_cur    = by_new

        S_now_xp     = _active_xp.abs(_active_xp.mean(_active_xp.exp(1j * theta_xp)))
        omega_max_xp = _active_xp.max(_active_xp.abs(omega_xp))
        mx_mean_xp   = _active_xp.mean(_active_xp.cos(theta_xp))
        my_mean_xp   = _active_xp.mean(_active_xp.sin(theta_xp))

        S_now     = float(S_now_xp)
        omega_max = float(omega_max_xp)
        mx_mean   = float(mx_mean_xp)
        my_mean   = float(my_mean_xp)

        _need_cpu_sync = (callback is not None) or (step % 20 == 0) or \
                         (frame_dir is not None and step % frame_every == 0)
        if _need_cpu_sync:
            theta_cur = _local_to_cpu(theta_xp).reshape(N, M)
            omega_cur = _local_to_cpu(omega_xp).reshape(N, M)

        if callback:
            callback(step, theta_cur.copy(), omega_cur.copy())

        if step % 20 == 0:
            hist.append((theta_cur.copy(), omega_cur.copy()))

        M_proj   = mx_mean * cos_phi + my_mean * sin_phi
        B_scalar = bx_cur * cos_phi + by_cur * sin_phi
        entry    = (t_now, B_scalar, M_proj, S_now)
        field_log.append(entry)
        if hysteresis_log is not None:
            hysteresis_log.append(entry)
        if sine_log is not None:
            sine_log.append(entry)

        if show_progress:
            _perf_now = _time_module.perf_counter()
            if (_perf_now - _perf_last_print_time) >= 0.15 or step == n_steps:
                _steps_since = step - _perf_last_print_step
                _dt_wall     = _perf_now - _perf_last_print_time
                _steps_per_s = _steps_since / _dt_wall if _dt_wall > 0 else 0.0
                _backend_tag = "GPU" if _active_gpu else "CPU"
                _frac_done   = step / n_steps if n_steps > 0 else 1.0
                B_status = bx_cur * cos_phi + by_cur * sin_phi
                _print_progress_bar(
                    _frac_done,
                    step=step, n_steps=n_steps, steps_per_s=_steps_per_s,
                    t_now=t_now, B_status=B_status, S_now=S_now, omega_max=omega_max,
                    backend=_backend_tag)
                _perf_last_print_step = step
                _perf_last_print_time = _perf_now

        _status_now = _time_module.perf_counter()
        if (_status_now - _last_status_print_time) >= 2.0 or step == n_steps:
            B_status = bx_cur * cos_phi + by_cur * sin_phi
            if show_progress and _progress_ansi_ok():
                _backend_tag = "GPU" if _active_gpu else "CPU"
                _frac_done   = step / n_steps if n_steps > 0 else 1.0
                _steps_since = step - _perf_last_print_step
                _dt_wall     = _status_now - _perf_last_print_time
                _steps_per_s = _steps_since / _dt_wall if _dt_wall > 0 else 0.0
                _print_progress_bar(
                    _frac_done,
                    step=step, n_steps=n_steps, steps_per_s=_steps_per_s,
                    t_now=t_now, B_status=B_status, S_now=S_now, omega_max=omega_max,
                    backend=_backend_tag)
            else:
                _status_text = (f"  [step {step}/{n_steps}]  t={t_now:.4f}s  "
                                f"B={B_status*1e3:+.4f}mT  S={S_now:.4f}  "
                                f"w_max={omega_max:.3f}rad/s")
                if show_progress:
                    _print_progress_bar_finish()
                _print(_status_text)
            _last_status_print_time = _status_now

        if frame_dir is not None and step % frame_every == 0:
            B_signed = bx_cur * cos_phi + by_cur * sin_phi
            _save_frame(step, theta_cur, omega_cur,
                        B_ext_inst=B_signed, phi_ext_inst=phi_ext_deg)

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

        if field_mode == 'pulse':
            if _pulse_phase == 'delay' and t_now >= field_delay:
                _pulse_phase = 'field_on'
                _finish_progress_bar_if_shown()
                _print(f"\n  Pulso: espera inicial concluída em t={t_now:.4f}s  campo LIGADO")
                tau_xp = _torques_xp(theta_xp, *field_at(t_now))
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

        if (field_mode == 'pulse' and _pulse_phase == 'field_on'
                and t_pulse is None):
            if S_now >= 0.99:
                _S99_count += 1
            else:
                _S99_count = 0

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
                _pulse_relaxing[0] = True
                _pulse_phase = 'relaxing'
                _finish_progress_bar_if_shown()
                if _trigger_A:
                    _print(f"\n  Pulso: S>=0.99 em t={t_now:.4f}s  campo zerado, relaxando")
                else:
                    _print(f"\n  Pulso: S estabilizou (S={S_now:.4f}) em t={t_now:.4f}s  campo zerado, relaxando")
                tau_xp = _torques_xp(theta_xp, 0.0, 0.0)
                if frame_dir is not None:
                    theta_cur = _local_to_cpu(theta_xp).reshape(N, M)
                    omega_cur = _local_to_cpu(omega_xp).reshape(N, M)
                    _save_frame(step, theta_cur, omega_cur,
                                stop_label="campo zerado - relaxando",
                                B_ext_inst=0.0)

        # ── condition 2: lattice in equilibrium ─────────────────────────────────
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
                _print(f"\n  Lattice at rest at t={t_now:.4f}s  (step {step}/{n_steps})  S={S_now:.4f}")
                _stop_reason = f"lattice at rest at t = {t_now:.4f} s"
                _stop_flag.set()
                if frame_dir is not None:
                    theta_cur = _local_to_cpu(theta_xp).reshape(N, M)
                    omega_cur = _local_to_cpu(omega_xp).reshape(N, M)
                    _save_frame(step, theta_cur, omega_cur,
                                stop_label="● lattice at rest")
                break

    _finish_progress_bar_if_shown()

    theta_cur = _local_to_cpu(theta_xp).reshape(N, M)
    omega_cur = _local_to_cpu(omega_xp).reshape(N, M)

    if _active_gpu:
        cp.cuda.Stream.null.synchronize()
    _perf_t_end       = _time_module.perf_counter()
    _perf_total_s     = _perf_t_end - _perf_t_start
    _perf_steps_done  = step
    _perf_steps_per_s = _perf_steps_done / _perf_total_s if _perf_total_s > 0 else 0.0
    _perf_ms_per_step = (_perf_total_s / _perf_steps_done * 1000.0) if _perf_steps_done > 0 else 0.0
    _perf_backend_tag = "GPU" if _active_gpu else "CPU"
    _print()
    _print(f"  ── performance ({_perf_backend_tag}) ──────────────────────────")
    _print(f"  Integration time    : {_perf_total_s:.3f} s  (wall-clock)")
    _print(f"  Steps executed      : {_perf_steps_done}")
    _print(f"  Throughput          : {_perf_steps_per_s:.1f} steps/s  "
           f"({_perf_ms_per_step:.4f} ms/step)")
    _print(f"  Needles in lattice  : {N*M}  (K = N x M)")

    _stop_flag.set()
    _restore_terminal()

    if hysteresis_log:
        import csv
        csv_path = "hysteresis_loop.csv"
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['t_s', 'B_T', 'M_proj', 'S'])
            w.writerows(hysteresis_log)
        _print(f"  Hysteresis data saved: {csv_path}")
        if make_images:
            _plot_hysteresis(hysteresis_log, hyst_autoscale=hyst_autoscale)

    if sine_log:
        import csv
        csv_path = "sine_field.csv"
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['t_s', 'B_T', 'M_proj', 'S'])
            w.writerows(sine_log)
        _print(f"  Sinusoidal field data saved: {csv_path}")
        if make_images:
            _plot_sine(sine_log, field_freq)

    return theta_cur, omega_cur, hist, n_frames, dt, _stop_reason, field_log


def next_available_path(path):
    import os, re
    base, ext = os.path.splitext(path)
    prefix = re.sub(r'\d{4}$', '', base)
    n = 0
    while True:
        candidate = f"{prefix}{n:04d}{ext}"
        if not os.path.exists(candidate):
            return candidate
        n += 1


def render_video(frame_dir, output_path, fps=24, crf=20, use_gpu=False):
    import subprocess, shutil, os

    if not shutil.which('ffmpeg'):
        _print("WARNING: ffmpeg not found.")
        return False

    input_pattern = os.path.join(frame_dir, "frame_%05d.png")
    vf = 'pad=ceil(iw/2)*2:ceil(ih/2)*2'

    strategies = []
    if use_gpu and _GPU_AVAILABLE:
        strategies.append((
            "H.264 via NVENC (GPU)",
            ['-c:v', 'h264_nvenc', '-preset', 'p4', '-cq', str(crf),
             '-pix_fmt', 'yuv420p']))

    strategies += [
        ("H.264 (libx264) — CPU",
         ['-c:v', 'libx264', '-preset', 'fast', '-crf', str(crf),
          '-pix_fmt', 'yuv420p']),
        ("MPEG-4 (mpeg4) — fallback for old ffmpeg",
         ['-c:v', 'mpeg4', '-q:v', '5',
          '-pix_fmt', 'yuv420p']),
        ("default ffmpeg codec",
         ['-pix_fmt', 'yuv420p']),
    ]

    _print(f"\nAssembling MP4 video: {output_path}")
    _print(f"  Source : {frame_dir}/frame_*.png   FPS={fps}   CRF={crf}")

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
            _print(f"  Saved : {output_path}  ({size_mb:.1f} MB)")
            return True
        else:
            last_err = result.stderr.strip().splitlines()
            short_err = last_err[-1] if last_err else "(no message)"
            _print(f"  [{desc}] failed: {short_err}")

    _print(f"\nERROR: ffmpeg failed to assemble the video.")
    return False


# ══════════════════════════════════════════════════════════════════════════════
# 5. GENERATION OF THE NEEDLE GRID
# ══════════════════════════════════════════════════════════════════════════════

def make_grid(N=8, M=8, geometry='square', noise=1.5, R=0.025):
    s3 = np.sqrt(3.0)

    if geometry == 'square':
        d = 2.0 * R
        xs = np.zeros((N, M))
        ys = np.zeros((N, M))
        for i in range(N):
            for j in range(M):
                xs[i, j] = j * d
                ys[i, j] = i * d
        thetas = noise * np.random.randn(N, M)
        Lx, Ly = M * d, N * d
        return xs, ys, thetas, d, Lx, Ly

    elif geometry == 'triangular':
        d      = 2.0 * R
        dx_col = d
        dy_row = R * s3
        offset = R
        xs = np.zeros((N, M))
        ys = np.zeros((N, M))
        for i in range(N):
            x_off = offset * (i % 2)
            for j in range(M):
                xs[i, j] = j * dx_col + x_off
                ys[i, j] = i * dy_row
        thetas = noise * np.random.randn(N, M)
        Lx, Ly = M * dx_col, N * dy_row
        return xs, ys, thetas, d, Lx, Ly

    elif geometry == 'honeycomb':
        dy  = R * np.sqrt(3.0)
        d   = 2.0 * R

        W = (M - 1) * 2.0 * R
        H = (N - 1) * dy

        N_rows = (N + 4) * 2
        x_start = -2.0 * 2 * R
        y_start = -2.0 * dy

        xs_list, ys_list = [], []
        for row in range(N_rows):
            y = y_start + row * dy
            if row % 2 == 0:
                x = x_start
                while x <= W + 2 * 2 * R:
                    xs_list.append(x)
                    ys_list.append(y)
                    x += 2.0 * R
            else:
                x = x_start + R
                while x <= W + 2 * 2 * R:
                    xs_list.append(x)
                    ys_list.append(y)
                    x += 4.0 * R

        all_x = np.array(xs_list)
        all_y = np.array(ys_list)

        margin = R * 0.99
        mask = ((all_x >= -margin) & (all_x <= W + margin) &
                (all_y >= -margin) & (all_y <= H + margin))
        clipped_x = all_x[mask]
        clipped_y = all_y[mask]

        n_pts = len(clipped_x)
        xs     = clipped_x.reshape(n_pts, 1)
        ys     = clipped_y.reshape(n_pts, 1)
        thetas = noise * np.random.randn(n_pts, 1)
        Lx, Ly = W, H
        return xs, ys, thetas, d, Lx, Ly

    else:
        raise ValueError(f"Geometria desconhecida: '{geometry}'")


# ══════════════════════════════════════════════════════════════════════════════
# 5b. IDENTIFICATION OF MAGNETIC DOMAINS
# ══════════════════════════════════════════════════════════════════════════════

def label_magnetic_domains(thetas, tol_deg=15.0):
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
        d = a - b
        return np.abs(np.arctan2(np.sin(d), np.cos(d)))

    idx_grid = np.arange(K).reshape(N, M)

    if N > 1:
        d_vert = _angdiff(thetas[:-1, :], thetas[1:, :])
        mask_vert = d_vert <= tol_rad
        pairs_a = [idx_grid[:-1, :][mask_vert]]
        pairs_b = [idx_grid[1:, :][mask_vert]]
    else:
        pairs_a, pairs_b = [], []

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
    if n_domains <= 1:
        return np.array([[0.55, 0.55, 0.85, 1.0]])
    hues = np.linspace(0.0, 1.0, n_domains, endpoint=False)
    hues = (hues * 0.618034 + 0.15) % 1.0
    return cm.hsv(hues)


# ══════════════════════════════════════════════════════════════════════════════
# 6. DESENHO OF A NEEDLE OF COMPASS
# ══════════════════════════════════════════════════════════════════════════════

def draw_compass(ax, x, y, theta, length=0.42, width=0.10,
                 color_n='#FFFFFF', color_s='#2E6DB4',
                 edge='#1a1a1a', zorder=4):
    half   = length / 2.0
    half_w = width  / 2.0

    pts_local = np.array([
        [ half,     0.0    ],
        [ 0.0,      half_w ],
        [-half,     0.0    ],
        [ 0.0,     -half_w ],
    ])

    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s],
                  [s,  c]])
    pts = (R @ pts_local.T).T + np.array([x, y])

    north = plt.Polygon(
        [pts[0], pts[1], [x, y], pts[3]],
        closed=True, facecolor=color_n,
        edgecolor=edge, linewidth=0.5, zorder=zorder)

    south = plt.Polygon(
        [pts[2], pts[1], [x, y], pts[3]],
        closed=True, facecolor=color_s,
        edgecolor=edge, linewidth=0.5, zorder=zorder)

    ax.add_patch(south)
    ax.add_patch(north)

    ax.plot(x, y, 'o', ms=2.0, color='#555555',
            markeredgecolor='#222222', markeredgewidth=0.4,
            zorder=zorder + 1)


def draw_compass_batch(ax, xs, ys, thetas, length=0.42, width=0.10,
                       color_n='#FFFFFF', color_s='#2E6DB4',
                       edge='#1a1a1a', zorder=4):
    half   = length / 2.0
    half_w = width  / 2.0
    x_flat = np.asarray(xs).ravel()
    y_flat = np.asarray(ys).ravel()
    th_flat = np.asarray(thetas).ravel()
    K = x_flat.shape[0]

    pts_local = np.array([
        [ half,     0.0    ],
        [ 0.0,      half_w ],
        [-half,     0.0    ],
        [ 0.0,     -half_w ],
    ])

    c = np.cos(th_flat)
    s = np.sin(th_flat)
    R = np.empty((K, 2, 2))
    R[:, 0, 0] = c
    R[:, 0, 1] = -s
    R[:, 1, 0] = s
    R[:, 1, 1] = c

    pts_rot = np.einsum('kab,vb->kva', R, pts_local)
    pts_rot[:, :, 0] += x_flat[:, None]
    pts_rot[:, :, 1] += y_flat[:, None]

    centers = np.stack([x_flat, y_flat], axis=1)

    north_patches = []
    south_patches = []
    for k in range(K):
        p0, p1, p2, p3 = pts_rot[k]
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

    scatter_pins = ax.scatter(
        x_flat, y_flat, s=4.0, c='#555555',
        edgecolors='#222222', linewidths=0.4,
        zorder=zorder + 1)

    return pc_north, pc_south, scatter_pins


# ══════════════════════════════════════════════════════════════════════════════
# 7. LATTICE-STATE VISUALIZATION
# ══════════════════════════════════════════════════════════════════════════════

def _fmt_B(B):
    if B == 0:      return "0 T"
    if B >= 0.1:    return f"{B:.4f} T"
    if B >= 1e-4:   return f"{B*1e3:.4f} mT"
    return              f"{B*1e6:.2f} µT"


def draw_ext_field_on_lattice(ax, xs, ys, B_ext, phi_ext_deg,
                              needle_len, B_ext_max=None, B_signed=None,
                              color='#FFD700'):
    """Represents the external field in the image."""
    b_val = B_signed if B_signed is not None else B_ext
    
    B_max_ref = B_ext_max if (B_ext_max and B_ext_max > 0) else max(abs(b_val), 1e-30)
    
    if B_signed is not None:
        frac = (B_signed / B_max_ref + 1.0) / 2.0
    else:
        frac = min(abs(B_ext) / B_max_ref, 1.0)
    frac = max(0.0, min(frac, 1.0))
    field_on  = (abs(b_val) > 1e-30)

    # ── 2. panel field (ALWAYS visible) — corner upper right ──────
    xlim  = ax.get_xlim()
    ylim  = ax.get_ylim()
    xspan = xlim[1] - xlim[0]
    yspan = ylim[1] - ylim[0]

    bar_w = 0.20 * xspan
    bar_h = 0.018 * yspan
    arrow_zone_w = 0.07 * xspan

    box_w = bar_w + arrow_zone_w
    px = xlim[1] - 0.02 * xspan - box_w
    py_top = ylim[1] - 0.035 * yspan

    pad = needle_len * 0.25
    box_h = 0.140 * yspan + pad

    # semi-transparent background box
    ax.add_patch(plt.Rectangle(
        (px - pad * 0.3, py_top - box_h),
        box_w + pad, box_h,
        facecolor='#080818', edgecolor=color,
        linewidth=0.8, alpha=0.82, zorder=19,
        transform=ax.transData))

    # ── text: intensity and direction ────────────
    if abs(b_val) > 1e-12:
        sign_str = "+" if b_val >= 0 else "-"
        b_str = sign_str + _fmt_B(abs(b_val))
    else:
        b_str = "0 T  (desligado)" if not field_on else "0 T"
    phi_str = f"dir: {phi_ext_deg:.1f} graus"

    text_color = color if field_on else '#777777'
    
    # Label "B ext" (left)
    ax.text(px, py_top - 0.022 * yspan, "B ext", color=color, fontsize=7,
            alpha=0.8, fontfamily='monospace', ha='left', va='top', zorder=23,
            transform=ax.transData)
            
    # Field value (second line)
    ax.text(px, py_top - 0.048 * yspan, b_str, color=text_color,
            fontsize=9, fontweight='bold', fontfamily='monospace',
            ha='left', va='top', zorder=23, transform=ax.transData)
            
    # Direction (third line)
    ax.text(px, py_top - 0.072 * yspan, phi_str, color=text_color,
            fontsize=7, fontfamily='monospace', ha='left', va='top', zorder=23,
            transform=ax.transData)

    # ── horizontal intensity bar ───
    bar_y = py_top - 0.115 * yspan
    ax.add_patch(plt.Rectangle(
        (px, bar_y), bar_w, bar_h,
        facecolor='#252540', edgecolor='none',
        zorder=20, transform=ax.transData))
    if frac > 0.01:
        ax.add_patch(plt.Rectangle(
            (px, bar_y), bar_w * frac, bar_h,
            facecolor=color if field_on else '#555555',
            edgecolor='none', alpha=0.85, zorder=21,
            transform=ax.transData))
    ax.add_patch(plt.Rectangle(
        (px, bar_y), bar_w, bar_h,
        facecolor='none', edgecolor='#5A5A8A',
        linewidth=0.7, zorder=22, transform=ax.transData))

    # ── direction dial (circle + arrow): centered next to the text lines ──
    arrow_cx = px + box_w - pad - 0.01 * xspan
    arrow_cy = py_top - 0.048 * yspan
    r_dial   = needle_len * 0.55 * 1.30
    
    dial = plt.Circle((arrow_cx, arrow_cy), r_dial, facecolor='#080818',
                      edgecolor='#5A5A8A', lw=0.7, zorder=20)
    ax.add_patch(dial)
    
    f = b_val / B_max_ref if B_max_ref > 0 else 0.0
    f = max(-1.0, min(f, 1.0))
    
    if field_on and abs(f) > 1e-5:
        phi_arrow = (phi_ext_deg + 180.0) % 360.0 if f < 0 else phi_ext_deg
        arrow_len = abs(f) * (2.0 * r_dial)
        
        dx = arrow_len * np.cos(np.deg2rad(phi_arrow))
        dy = arrow_len * np.sin(np.deg2rad(phi_arrow))
        
        x_tail = arrow_cx - dx / 2.0
        y_tail = arrow_cy - dy / 2.0
        x_head = arrow_cx + dx / 2.0
        y_head = arrow_cy + dy / 2.0
        arrow_color = color
        ax.annotate(
            '', xy=(x_head, y_head), xytext=(x_tail, y_tail),
            arrowprops=dict(
                arrowstyle='->', color=arrow_color, lw=1.8,
                mutation_scale=10),
            zorder=21)
    else:
        ax.plot(arrow_cx, arrow_cy, 'o', ms=2, color='#555555', zorder=21)

    # ── 1. arrows at the sites ─────────────────
    if field_on:
        phi_inst = np.deg2rad(phi_ext_deg + 180.0 if b_val < 0 else phi_ext_deg)
        cos_phi  = np.cos(phi_inst)
        sin_phi  = np.sin(phi_inst)
        alen_sites = needle_len * 0.35
        dx_a = alen_sites * cos_phi
        dy_a = alen_sites * sin_phi
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
    N, M = thetas.shape
    _r = r_halo if r_halo is not None else needle_len * 0.58
    _n_domains_found = None
    _n_significant   = None

    if halo_mode == 'domains':
        domain_labels, _n_domains_found = label_magnetic_domains(
            thetas, tol_deg=domain_tol_deg)
        K_total = N * M
        _unique_labels, _inverse, _counts = np.unique(
            domain_labels, return_inverse=True, return_counts=True)
        _inverse = _inverse.ravel()
        _min_domain_size = max(3, int(0.02 * K_total))
        _is_significant = _counts >= _min_domain_size
        _n_significant = int(np.sum(_is_significant))

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
        align_sum   = np.zeros((N, M))
        align_count = np.zeros((N, M))
        for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
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
    """Generates a figure with the instantaneous state of the needle lattice."""
    if figsize_inches is None:
        x_span = xs.max() - xs.min() + 2 * needle_len * 2.0
        y_span = ys.max() - ys.min() + 2 * needle_len * 4.0
        aspect = x_span / y_span if y_span > 0 else 1.0
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

    if show_order:
        _n_domains_found, _n_significant = draw_halos_batch(
            ax, xs, ys, thetas, needle_len, r_halo=r_halo,
            halo_mode=halo_mode, domain_tol_deg=domain_tol_deg)

    draw_compass_batch(ax, xs, ys, thetas,
                       length=needle_len, width=needle_width)

    ax.set_aspect('equal')
    margin = needle_len * 1.6
    yspan_lattice = ys.max() - ys.min() + 2 * margin
    yspan = yspan_lattice / 0.8
    ax.set_xlim(xs.min() - margin, xs.max() + margin)
    ax.set_ylim(ys.min() - margin, ys.min() - margin + yspan)
    ax.set_title(title, color='#ECF0F1', fontsize=11, pad=10,
                 fontfamily='monospace')
    ax.tick_params(colors='#7F8C8D')
    for spine in ax.spines.values():
        spine.set_edgecolor('#2C3E50')

    draw_ext_field_on_lattice(ax, xs, ys, B_ext, phi_ext_deg, needle_len,
                              B_ext_max=B_ext_max, B_signed=B_signed)

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
    """Generates a frame-by-frame animation of the lattice time evolution."""
    fig, ax = plt.subplots(figsize=(7, 7), facecolor='#1A1A2E')
    ax.set_facecolor('#16213E')
    ax.set_aspect('equal')
    margin = needle_len * 1.6
    yspan_lattice = ys.max() - ys.min() + 2 * margin
    yspan = yspan_lattice / 0.8
    ax.set_xlim(xs.min() - margin, xs.max() + margin)
    ax.set_ylim(ys.min() - margin, ys.min() - margin + yspan)

    N, M = thetas_hist[0].shape

    def update(frame):
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
        _print(f"GIF saved at: {save_gif}")
    return fig, ani


# ══════════════════════════════════════════════════════════════════════════════
# 9. GLOBAL ORDER-PARAMETER PLOT
# ══════════════════════════════════════════════════════════════════════════════

def plot_order_parameter(thetas_hist, outpath, dt=None):
    """Plots the time evolution of the global magnetic order parameter S(t)."""
    order_params = [
        np.abs(np.mean(np.exp(1j * th)))
        for th in thetas_hist
    ]

    if dt is not None:
        time_ax = np.arange(len(order_params)) * 20 * dt
        xlabel  = "Time  t  [s]"
    else:
        time_ax = np.arange(len(order_params))
        xlabel  = "Snapshot (every 20 steps)"

    fig, ax = plt.subplots(figsize=(8, 3.5), facecolor='#1A1A2E')
    ax.set_facecolor('#0F3460')
    ax.plot(time_ax, order_params, color='#E94560', lw=2)
    ax.set_xlabel(xlabel, color='#BDC3C7')
    ax.set_ylabel(r"Order parameter $S = |\langle e^{i\theta}\rangle|$",
                  color='#BDC3C7')
    ax.set_title("Global Magnetic Order Parameter Evolution",
                 color='#ECF0F1', fontfamily='monospace')
    ax.tick_params(colors='#7F8C8D')
    ax.set_ylim(0, 1.05)
    ax.grid(True, color='#2C3E50', alpha=0.5)
    for spine in ax.spines.values():
        spine.set_edgecolor('#2C3E50')
    plt.tight_layout()
    plt.savefig(outpath, dpi=130, bbox_inches='tight', facecolor='#1A1A2E')
    _print(f"Order parameter saved at: {outpath}")


# ══════════════════════════════════════════════════════════════════════════════
# 10. INTERFACE OF LINE OF COMANDO (main)
# ══════════════════════════════════════════════════════════════════════════════

def main():
    """Command-line entry point of the simulation."""
    parser = argparse.ArgumentParser(
        description="Simulacao de rede de bussolas — campo dipolar 2D",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    grp = parser.add_argument_group("Geometria da rede")
    grp.add_argument('--N', type=int, default=8,
                     help='Numero de LINHAS de agulhas')
    grp.add_argument('--M', type=int, default=8,
                     help='Numero de COLUNAS de agulhas')
    grp.add_argument('--R', type=float, default=0.025,
                     help='Raio do circulo que envolve cada agulha [m]. '
                           'A distancia entre centers de vizinhos = 2R. '
                           'Padrao: 0.025 m = 2.5 cm')
    grp.add_argument('--needle_frac', type=float, default=0.80,
                     help='Comprimento da agulha como fracao do diametro 2R '
                           '(0.0 a 0.8). Padrao: 0.80  ->  agulha = 0.80 * 2R. '
                           'Valores fora do intervalo serao limitados (clamp).')
    grp.add_argument('--geometry',
                     choices=['square', 'triangular', 'honeycomb'],
                     default='square',
                     help='Tipo de rede')

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
                           'inercia. Padrao: {STEEL_DENSITY_DEFAULT:.0f} '
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
                      help='Modo do campo externo')
    grp2.add_argument('--hyst_spacing', choices=['linear', 'log'],
                      default='linear',
                      help='Espacamento da rampa de histerese')
    grp2.add_argument('--hyst_log_k', type=float, default=5.0,
                      help='Concentracao do espacamento log da histerese')
    grp2.add_argument('--field_delay', type=float, default=0.0,
                      help='Delta_t de espera [s]')
    grp2.add_argument('--t_pulse', type=float, default=None,
                      help='Duracao do pulso [s]')
    grp2.add_argument('--torque_tol', type=float, default=1e-3,
                      help='Tolerancia relativa do criterio de equilibrio por torque')
    grp2.add_argument('--field_freq', type=float, default=1.0,
                      help='Frequencia do campo senoidal [Hz].')
    grp2.add_argument('--dt_factor', type=float, default=0.05,
                      help='Fracao do periodo natural T0 usada como passo dt')
    grp2.add_argument('--noise', type=float, default=1.5,
                      help='Amplitude do ruido inicial nos angulos [rad]')
    grp2.add_argument('--seed', type=int, default=42,
                      help='Semente do gerador aleatorio')
    grp2.add_argument('--pbc', type=int, choices=[0, 1], default=0,
                      help='Condicoes periodicas de contorno (PBC)')
    grp2.add_argument('--pbc_images', type=int, default=1,
                      help='Numero de replicas periodicas somadas')
    grp2.add_argument('--gpu', type=int, choices=[0, 1], default=0,
                      help='Uso de GPU (CuPy)')
    grp2.add_argument('--progress_bar', type=int, choices=[0, 1], default=1,
                      help='Exibe a barra de progresso em tempo real')
    grp2.add_argument('--halo_mode', type=str, choices=['order', 'domains'],
                      default='order',
                      help='Modo de coloracao dos halos')
    grp2.add_argument('--domain_tol', type=float, default=15.0,
                      help='Tolerancia angular [graus] para dominios')

    grp3 = parser.add_argument_group("Campo externo uniforme (SI)")
    grp3.add_argument('--B_ext', type=float, default=0.0,
                      help='Intensidade do campo externo [T].')
    grp3.add_argument('--phi_ext', type=float, default=0.0,
                      help='Direcao do campo externo [graus]')
    grp3.add_argument('--ext_Bx', type=float, default=None,
                      help='Componente Bx do campo externo [T]')
    grp3.add_argument('--ext_By', type=float, default=None,
                      help='Componente By do campo externo [T]')

    grp4 = parser.add_argument_group("Saida")
    grp4.add_argument('--video', type=str, default=None,
                      metavar='NOME',
                      help='Gera video MP4 com este nome base.')
    grp4.add_argument('--frame_every', type=int, default=5,
                      help='Salva um frame a cada N passos.')
    grp4.add_argument('--make_images', type=int, default=1, choices=[0, 1],
                      help='1 (padrao) = gera as imagens PNG, 0 = desliga.')
    grp4.add_argument('--fps', type=int, default=24,
                      help='Quadros por segundo do video MP4.')
    grp4.add_argument('--dpi', type=int, default=120,
                      help='Resolucao dos frames em DPI.')
    grp4.add_argument('--keep_frames', action='store_true',
                      help='Mantém a pasta de PNGs intermediários')
    grp4.add_argument('--hyst_autoscale', type=int, choices=[0, 1], default=0,
                      help='Enables automatic scaling (autoscale) of the axes for the MxH hysteresis curve.')
    grp4.add_argument('--csv_order', choices=['t', 'B'], default='t',
                      help='Ordem das colunas no CSV exportado.')
    args = parser.parse_args()

    if args.needle_frac < 0.0 or args.needle_frac > 0.8:
        clamped = max(0.0, min(args.needle_frac, 0.8))
        _print(f"  Aviso: --needle_frac {args.needle_frac} fora do intervalo "
                f"[0, 0.8]; ajustado para {clamped}")
        args.needle_frac = clamped

    np.random.seed(args.seed)

    if args.ext_Bx is not None or args.ext_By is not None:
        Bext_x = args.ext_Bx if args.ext_Bx is not None else 0.0
        Bext_y = args.ext_By if args.ext_By is not None else 0.0
        B_ext_mag = np.sqrt(Bext_x**2 + Bext_y**2)
        phi_ext_deg = np.degrees(np.arctan2(Bext_y, Bext_x))
    else:
        phi_rad = np.deg2rad(args.phi_ext)
        Bext_x  = args.B_ext * np.cos(phi_rad)
        Bext_y  = args.B_ext * np.sin(phi_rad)
        B_ext_mag   = args.B_ext
        phi_ext_deg = args.phi_ext

    R          = args.R
    needle_len = args.needle_frac * 2.0 * R
    needle_width = needle_len * 0.22

    if args.inertia is None:
        inertia = compute_inertia_from_geometry(
            needle_len, needle_width, args.needle_thickness,
            density=args.steel_density)
        _inertia_auto = True
    else:
        inertia = args.inertia
        _inertia_auto = False

    if args.steel_Bsat is not None:
        _Ms_used = args.steel_Bsat / (4.0 * np.pi * 1e-7)
    else:
        _Ms_used = STEEL_MS_SATURATION_DEFAULT
    if args.moment is None:
        moment = compute_moment_from_geometry(
            needle_len, needle_width, args.needle_thickness, Ms=_Ms_used)
        _moment_auto = True
    else:
        moment = args.moment
        _moment_auto = False

    def fmt_field(B):
        if B == 0:      return "0 T"
        if B >= 0.1:    return f"{B:.4f} T"
        if B >= 1e-4:   return f"{B*1e3:.4f} mT"
        return              f"{B*1e6:.2f} µT"

    _print(f"\n{'═'*62}")
    _print(f"  Version      : V55")
    _print(f"  Developer    : J.P. Sinnecker")
    if args.gpu and _GPU_AVAILABLE:
        try:
            _gpu_name = cp.cuda.runtime.getDeviceProperties(0)['name'].decode()
        except Exception:
            _gpu_name = "GPU CUDA"
        _print(f"  Backend      : GPU ({_gpu_name}) via CuPy  [--gpu 1]")
    elif args.gpu and not _GPU_AVAILABLE:
        _print(f"  Backend      : CPU (NumPy)  [--gpu 1 requested, but GPU unavailable]")
        if _GPU_ERROR_MSG:
            if 'CUDA headers' in _GPU_ERROR_MSG or 'CUDA_PATH' in _GPU_ERROR_MSG:
                _print(f"  GPU unavailable: missing CUDA Toolkit headers.")
                _print(f"  Solution: pip install cupy-cuda12x[ctk]")
            else:
                _print(f"  GPU unavailable: {_GPU_ERROR_MSG[:70]}")
        else:
            _print(f"  GPU unavailable: cupy not installed")
    else:
        _gpu_hint = " (GPU detected, use --gpu 1 to activate)" if _GPU_AVAILABLE else ""
        _print(f"  Backend      : CPU (NumPy)  [--gpu 0]{_gpu_hint}")
    _print(f"  Array        : {args.geometry}  {args.N}x{args.M} needles")
    pbc_str = (f"on  (n_images={args.pbc_images}, sum over "
               f"{(2*args.pbc_images+1)**2} cells)") if args.pbc else "off"
    _print(f"  PBC          : {pbc_str}")
    _print(f"  Radius R     : {args.R*100:.2f} cm  (2R = {2*args.R*100:.2f} cm)")
    _print(f"  Needle R     : {args.needle_frac*100:.0f}% of 2R = {args.needle_frac*2*args.R*100:.2f} cm")
    if _moment_auto:
        _Bsat_used = _Ms_used * 4.0 * np.pi * 1e-7
        _print(f"  Mag Moment.  : {moment:.4g} A.m2  per needle")
    else:
        _print(f"  Mag Moment.  : {moment:.4g} A.m2  per needle [manual value]")
    if _inertia_auto:
        _mass_g = args.steel_density * needle_len * needle_width * args.needle_thickness * 1e3
        _print(f"  Inertia      : {inertia:.3e} kg.m2  per needle")
    else:
        _print(f"  Inertia      : {inertia:.3e} kg.m2  per needle [manual value]")
    _print(f"  Damping      : {args.damping:.3e} N.m.s/rad  (air)")
    _print(f"  External Field: {fmt_field(B_ext_mag)}  phi={phi_ext_deg:.1f} degrees")
    _print(f"  Components   : Bx={fmt_field(abs(Bext_x))}  By={fmt_field(abs(Bext_y))}")
    field_mode_str = args.field_mode
    _print(f"  Field Mode   : {field_mode_str}")
    _print(f"  Noise Seed   : {args.seed}")
    _print(f"{'═'*62}\n")

    xs, ys, thetas_init, nn_dist, Lx_period, Ly_period = make_grid(
        N=args.N, M=args.M,
        geometry=args.geometry,
        noise=args.noise,
        R=args.R,
    )

    r_halo     = R * 0.98
    cutoff     = nn_dist * 2.6

    if args.pbc:
        max_cutoff_pbc = min(Lx_period, Ly_period) / 2.0
        if cutoff > max_cutoff_pbc:
            _print(f"  PBC: cutoff reduzido de {cutoff*100:.2f}cm para "
                    f"{max_cutoff_pbc*100:.2f}cm (limite min(Lx,Ly)/2)")
            cutoff = max_cutoff_pbc

    B_ref = MU0_OVER_4PI * 2.0 * moment / nn_dist**3
    _print(f"  B_dipolar ref: {fmt_field(B_ref)}  (Between neighbours)")
    if B_ext_mag > 0:
        ratio = B_ext_mag / B_ref
        dom = "DOMINANT" if ratio > 1 else "weak"
        _print(f"  B_ext/B_ref  : {ratio:.3f}  (External Field - {dom})")
    _print()

    ext_kwargs = dict(B_ext=B_ext_mag, phi_ext_deg=phi_ext_deg,
                      r_halo=r_halo, halo_mode=args.halo_mode,
                      domain_tol_deg=args.domain_tol)

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

    if args.make_images:
        fig0, _ = plot_state(thetas_init, xs, ys,
                             title="Estado inicial (aleatório)",
                             needle_len=needle_len, needle_width=needle_width,
                             **ext_kwargs)
        plt.tight_layout()
        plt.savefig("compass_initial.png",
                    dpi=args.dpi, bbox_inches='tight', facecolor='#1A1A2E')
        plt.close(fig0)
        _print("Initial State Saved.")

    frame_dir   = None
    final_video = None
    if args.video and not args.make_images:
        _print(f"  Warning: --video '{args.video}' was requested along with "
               f"--make_images 0. Ignoring --video.")
        args.video = None

    if args.video:
        import os
        if not os.path.splitext(args.video)[1]:
            args.video = args.video + ".mp4"
        final_video = next_available_path(args.video)
        _print(f"  Video filename: '{final_video}'")
        base      = os.path.splitext(final_video)[0]
        frame_dir = base + "_frames"
        _print(f"Integrating and saving frames at '{frame_dir}/'...")
    else:
        _print("Integrating inertial dynamics...")

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
        hyst_autoscale=bool(args.hyst_autoscale),
    )
    frames_str = f"  ({n_frames} frames saved)" if n_frames else ""
    _print(f"Integration completed - {stop_reason}{frames_str}")

    import csv as _csv, os as _os
    if args.video and final_video:
        csv_path = _os.path.splitext(final_video)[0] + ".csv"
    else:
        csv_path = "compass_field_log.csv"

    if args.csv_order == 'B':
        header = ['B_applied_T', 'M_proj', 'S', 't_s']
        rows   = [(B, Mp, S, t) for (t, B, Mp, S) in field_log]
    else:
        header = ['t_s', 'B_applied_T', 'M_proj', 'S']
        rows   = field_log

    with open(csv_path, 'w', newline='') as f:
        w = _csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    _print(f"  CSV saved: {csv_path}  ({len(field_log)} points)")

    if args.make_images:
        title_eq = (f"Equilibrium — {args.geometry} {args.N}×{args.M}"
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
        _print("Equilibrium state saved.")

        fig2, axes = plt.subplots(1, 2, figsize=(15, 7), facecolor='#1A1A2E')
        margin = needle_len * 1.6
        yspan_lattice = ys.max() - ys.min() + 2 * margin
        yspan = yspan_lattice / 0.8
        for ax_ in axes:
            ax_.set_facecolor('#16213E')
            ax_.set_aspect('equal')
            ax_.set_xlim(xs.min() - margin, xs.max() + margin)
            ax_.set_ylim(ys.min() - margin, ys.min() - margin + yspan)
            ax_.tick_params(left=False, bottom=False,
                            labelleft=False, labelbottom=False)
            for sp in ax_.spines.values():
                sp.set_edgecolor('#2C3E50')

        axes[0].set_title("Initial State", color='#ECF0F1',
                          fontsize=12, fontfamily='monospace')
        axes[1].set_title("Dipolar Equilibrium", color='#ECF0F1',
                          fontsize=12, fontfamily='monospace')

        N_g, M_g = thetas_init.shape
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
                f"Initial State — {_ns0} main domains "
                f"(out of {_nd0})", color='#ECF0F1', fontsize=12,
                fontfamily='monospace')
            axes[1].set_title(
                f"Dipolar Equilibrium — {_ns1} main domains "
                f"(out of {_nd1})", color='#ECF0F1', fontsize=12,
                fontfamily='monospace')

        draw_compass_batch(axes[0], xs, ys, thetas_init,
                           length=needle_len, width=needle_width)
        draw_compass_batch(axes[1], xs, ys, thetas_eq,
                           length=needle_len, width=needle_width)

        draw_ext_field_on_lattice(axes[0], xs, ys, B_ext_mag, phi_ext_deg, needle_len)
        draw_ext_field_on_lattice(axes[1], xs, ys, B_ext_mag, phi_ext_deg, needle_len)

        bfield_str = (f"  |  B_ext={fmt_field(B_ext_mag)} @ {phi_ext_deg:.0f}°"
                      if B_ext_mag > 0 else "")
        fig2.suptitle(
            f"Array {args.geometry} {args.N}×{args.M}"
            f"  |  R={args.R*100:.1f} cm  2R={2*args.R*100:.1f} cm"
            f"{bfield_str}  |  2D dipolar interaction",
            color='#BDC3C7', fontsize=11, fontfamily='monospace')
        plt.tight_layout()
        plt.savefig("compass_comparison.png",
                    dpi=args.dpi, bbox_inches='tight', facecolor='#1A1A2E')
        plt.close(fig2)
        _print("Comparison saved.")

        plot_order_parameter([th for th, _ in thetas_hist],
                             "compass_order_param.png",
                             dt=sim_dt)

    if final_video and frame_dir and n_frames > 0:
        import shutil
        ok = render_video(frame_dir, final_video, fps=args.fps, use_gpu=bool(args.gpu))
        if ok and not args.keep_frames:
            shutil.rmtree(frame_dir)
            _print(f"  Frames folder removed: {frame_dir}/")
        elif not ok:
            _print(f"  PNG frames kept at: {frame_dir}/")

    _print("\nCompleted.")


if __name__ == '__main__':
    main()
