"""
============================================================
compass_sim.py — Compass needle lattice simulation
============================================================

Models a 2D grid of classical magnetic dipoles (compass
needles) that interact via the field each generates on its neighbors.
The dynamics are inertial (Newton's 2nd law of rotation) without friction on the pin,
with viscous air damping. Velocity-Verlet integrator.
All physical quantities in SI units.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
COMMAND LINE PARAMETERS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

┌─ LATTICE GEOMETRY ─────────────────────────────────────────┐
│                                                            │
│  --geometry  square | triangular | honeycomb               │
│              Lattice type. Default: square                 │
│              · square     : rectangular grid               │
│              · triangular : equilateral triangular grid    │
│              · honeycomb  : honeycomb (hexagonal holes)    │
│                                                            │
│  --N         int   Number of ROWS of needles. Default: 8   │
│  --M         int   Number of COLUMNS of needles. Default: 8│
│                                                            │
│  --R         float Radius of each needle's circle [m].     │
│                    Distance between neighbors = 2R.        │
│                    Default: 0.025 m (2.5 cm)               │
│                                                            │
│  --needle_frac float Length of the needle as a fraction    │
│                    of the diameter 2R (0.0–1.0).           │
│                    Default: 0.80  → needle = 80% of 2R     │
└────────────────────────────────────────────────────────────┘

┌─ PHYSICS AND SIMULATION ───────────────────────────────────┐
│                                                            │
│  --moment    float Magnetic moment of each needle [A·m²]   │
│                    Default: 0.1  (desk compass, ~5 cm)     │
│                    Ref: pocket ≈ 0.01 | nautical ≈ 1.0     │
│                                                            │
│  --inertia   float Moment of inertia [kg.m2]. If omitted,  │
│                    calculated automatically from geometry  │
│                    (steel blade: R, needle_frac,           │
│                    --needle_thickness, --steel_density)    │
│                                                            │
│  --damping   float Viscous air damping [N·m·s/rad]         │
│                    Controls the quality factor Q:          │
│                    Q = omega_0·I/b (High Q = more oscill.) │
│                    Default: 5e-8 (Q≈25, realistic compass) │
│                    For smooth B_ext=0.1T: use 8e-6 (Q≈4)   │
│                                                            │
│  --t_sim     float Total physical simulation time [s].     │
│                    Sum of all integrated dt steps          │
│                    — equivalent to what a real stopwatch   │
│                    would show watching the needles move.   │
│                    The video displays this physical time.  │
│                    Sim. stops early if S=1.00 or at rest.  │
│                    Default: 2.0 s                          │
│                                                            │
│  --dt_factor float Fraction of the natural period T₀ used  │
│                    as integration step (0.02–0.10).        │
│                    Smaller = more accurate, slower.        │
│                    Default: 0.05                           │
│                                                            │
│  --noise     float Initial noise amplitude [rad].          │
│                    0 = all point to +x                     │
│                    π ≈ 3.14 = completely random orient.    │
│                    Default: 1.5                            │
│                                                            │
│  --seed      int   Random generator seed.                  │
│                    Ensures reproducibility.                │
│                    Default: 42                             │
└────────────────────────────────────────────────────────────┘

┌─ UNIFORM EXTERNAL FIELD (SI units) ────────────────────────┐
│  Two ways to specify — do not use both together.           │
│                                                            │
│  Way A — intensity + angle (recommended):                  │
│  --B_ext     float Field intensity [T].                    │
│                    0.0      = no field (default)           │
│                    50e-6    = Earth's field (≈50 µT)       │
│                    1e-3     = fridge magnet at 5 cm        │
│                    0.1      = strong field (aligns all)    │
│                                                            │
│  --phi_ext   float Field direction [degrees].              │
│                    0   = right (+x)  ← default             │
│                    90  = up (+y)                           │
│                    180 = left (−x)                         │
│                    270 = down (−y)                         │
│                    Counter-clockwise. Default: 0.0         │
│                                                            │
│  Way B — cartesian components (overrides A):               │
│  --ext_Bx    float x component of the field [T]            │
│  --ext_By    float y component of the field [T]            │
└────────────────────────────────────────────────────────────┘

┌─ OUTPUT ───────────────────────────────────────────────────┐
│  PNG files always generated in the current directory:      │
│    compass_initial.png      initial state                  │
│    compass_equilibrium.png  final state                    │
│    compass_comparison.png   side-by-side comparison        │
│    compass_order_param.png  order parameter S(t)           │
│                                                            │
│  --video     str   Path of the MP4 video to generate.      │
│                    Requires ffmpeg installed.              │
│                    If the file already exists, saves as    │
│                    name0001.mp4, name0002.mp4, etc.        │
│                    Ex: --video simulation.mp4              │
│                                                            │
│  --frame_every int Saves a frame every N steps.            │
│                    Smaller = smoother video, slower.       │
│                    Default: 5                              │
│                                                            │
│  --fps       int   Frames per second of the MP4 video.     │
│                    Default: 24                             │
│                                                            │
│  --keep_frames     If present, keeps the intermediate PNG  │
│                    folder after generating the MP4.        │
└────────────────────────────────────────────────────────────┘

┌─ CONTROLS DURING SIMULATION ───────────────────────────────┐
│  The simulation automatically stops when:                  │
│    · S = 1.00  (all needles aligned)                       │
│    · Lattice at rest (ω_max → 0)                           │
│    · Time t_sim reached                                    │
│    · Ctrl+I (Tab) pressed in the terminal                  │
│      → stops and saves the video immediately               │
└────────────────────────────────────────────────────────────┘

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
EXAMPLES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  # Defaults — 8×8 square lattice without field
  python compass_sim.py

  # Honeycomb with Earth's field
  python compass_sim.py --geometry honeycomb --N 10 --M 10 --B_ext 50e-6

  # Triangular with 0.1 T field at 45°, smooth movement
  python compass_sim.py --geometry triangular --N 10 --M 10 \
      --B_ext 0.1 --phi_ext 45 --damping 8e-6 --t_sim 2.0

  # Video with larger needles and many oscillations
  python compass_sim.py --R 0.03 --needle_frac 0.85 --damping 1e-9 \
      --t_sim 5.0 --frame_every 2 --fps 30 --video sim.mp4

  # Field via cartesian components
  python compass_sim.py --ext_Bx 0.05 --ext_By -0.05

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Dependencies: numpy, matplotlib, ffmpeg (for video)
Optional (GPU acceleration): cupy-cuda12x (or version compatible with your CUDA)
  pip install cupy-cuda12x[ctk]
  The [ctk] suffix also installs the CUDA Toolkit headers, necessary
  to compile kernels at runtime (without them, falls back to CPU
  with error "Failed to find CUDA headers").
  Automatically detected; if missing or incomplete, uses CPU/NumPy.
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
[END_MARKER]
