# User Guide: Compass-Needle Lattice Simulator (`compass_sim`)

**Compiled Version:** July 8, 2026 at 10:42 (Markdown updated; PDF compiled dynamically)

Welcome to the **Compass-Needle Lattice Simulator** (`compass_sim`), a Python-based physics simulation engine that models the dynamic behavior of 2D arrays of interacting classical magnetic dipoles (compass needles). 

The simulation models real Newtonian rotational dynamics in SI units, integrating the equation of motion for each needle using the **Velocity-Verlet** algorithm. It is designed to run efficiently on both CPUs (via NumPy) and GPUs (via CuPy).

---

## Table of Contents
1. [Physics & Theoretical Model](#1-physics--theoretical-model)
2. [Lattice Geometry Options](#2-lattice-geometry-options)
3. [External Field Modes (`--field_mode`)](#3-external-field-modes)
4. [Demagnetization Protocols (`--demag`)](#4-demagnetization-protocols)
5. [System Requirements & Acceleration](#5-system-requirements--acceleration)
6. [Command Line Parameters Reference](#6-command-line-parameters-reference)
7. [Output Files & Data Formats](#7-output-files--data-formats)
8. [Step-by-Step Practical Examples](#8-step-by-step-practical-examples)
9. [Running Sweep Campaigns](#9-running-sweep-campaigns)
10. [Summary Reference Table of Parameters](#10-summary-reference-table-of-parameters)

---

## 1. Physics & Theoretical Model

At the core of the simulator is the rotational equation of motion for each needle $i$:

$$I \ddot{\theta}_i = \tau_i^{\text{dip}} + \tau_i^{\text{ext}} - b \dot{\theta}_i + \tau_i^{\text{noise}}$$

Where:
* **$\theta_i$**: Orientation angle of the $i$-th compass needle in the 2D plane ($xy$-plane).
* **$I$**: Moment of inertia of each needle ($kg \cdot m^2$).
* **$\tau_i^{\text{dip}}$**: Rotational torque due to dipole-dipole interactions with all other needles.
* **$\tau_i^{\text{ext}}$**: Rotational torque due to the applied uniform external magnetic field ($B_{\text{ext}}$).
* **$b$**: Viscous air damping coefficient ($N \cdot m \cdot s/\text{rad}$).
* **$\tau_i^{\text{noise}}$**: Random thermal or ambient torque fluctuations (active during simulated annealing).

### Automatic Property Calibration
If not supplied explicitly via the CLI, the magnetic moment $m$ ($A \cdot m^2$) and moment of inertia $I$ ($kg \cdot m^2$) are calculated automatically based on the geometry of the needle, assuming it is a flat rectangular steel sheet:
* **Mass**: Computed from thickness, length ($2R \cdot \text{needle\_frac}$), width ($L / 10$), and steel density ($\rho \approx 7850\text{ kg/m}^3$).
* **Inertia ($I$)**: $\frac{1}{12} \cdot \text{mass} \cdot (\text{length}^2 + \text{width}^2)$.
* **Magnetic Moment ($m$)**: $M_s \cdot \text{volume}$, where $M_s = B_{\text{sat}} / \mu_0$ (saturation magnetization of steel, default $B_{\text{sat}} = 2.0\text{ T}$).

---

## 2. Lattice Geometry Options

The layout of the needles affects the dipolar interaction field profile. The simulator supports three distinct 2D geometry layouts, selected using `--geometry`:

### 1. Square Lattice (`--geometry square`)
* A rectangular grid of dimensions $N \times M$.
* Nearest-neighbor spacing is $2R$.
* Leads to a highly symmetric anisotropic energy landscape.

### 2. Triangular Lattice (`--geometry triangular`)
* An equilateral triangular grid.
* Nearest-neighbor spacing is $2R$ along all directions at $60^{\circ}$ intervals.
* Introduces geometric frustration due to the triangular layout.

### 3. Honeycomb Lattice (`--geometry honeycomb`)
* A hexagonal grid layout (like graphene) with hexagonal voids.
* Characterized by a lower coordination number (3 nearest neighbors instead of 4 or 6).

---

## 3. External Field Modes (`--field_mode`)

You can control how the external magnetic field $B_{\text{ext}}(t)$ changes over time using `--field_mode`:

| Mode | Description | Key Parameters |
| :--- | :--- | :--- |
| `static` | Applies a constant field in a fixed direction. | `--B_ext`, `--phi_ext` (or `--ext_Bx`, `--ext_By`) |
| `hysteresis` | Sweeps the field through a full hysteresis cycle: $0 \rightarrow +B_{\text{max}} \rightarrow 0 \rightarrow -B_{\text{max}} \rightarrow 0 \rightarrow +B_{\text{max}}$. Spacing can be linear or logarithmic. | `--B_ext` (defines $B_{\text{max}}$), `--phi_ext`, `--hyst_spacing`, `--hyst_log_k` |
| `sine` | A sinusoidal external field: $B_{\text{ext}}(t) = B_{\text{max}} \sin(2\pi f t)$. | `--B_ext` (defines $B_{\text{max}}$), `--field_freq` |
| `pulse` | Applies a pulse of field magnitude `--B_ext` starting after `--field_delay` for a duration of `--t_pulse`. | `--B_ext`, `--field_delay`, `--t_pulse` |
| `step_pos` | Instantly steps the field from $0$ to $+B_{\text{max}}$ at $t = 0$. | `--B_ext` |
| `step_neg` | Instantly steps the field from $0$ to $-B_{\text{max}}$ at $t = 0$. | `--B_ext` |
| `forc` | First-Order Reversal Curves (FORC) protocol. Ramps to positive saturation, reverses to a varying reversal field $B_r$, then sweeps back up to saturation, recording minor loops to map hysteretic components. | `--forc_n_curves`, `--forc_Br_min`, `--forc_rate`, `--forc_t_ramp_up`, `--forc_t_ramp_down` |

---

## 4. Demagnetization Protocols (`--demag`)

Demagnetization is used to initialize the lattice in a randomized state before starting the primary field simulation. Specify the protocol via `--demag`:

### 1. Biaxial Rotational Demagnetization (`--demag on` or `--demag rotational`)
* **Mechanism**: The external field vector rotates continuously in the $xy$-plane while its magnitude decays linearly to zero.
  $$B_x(t) = B_{\text{max}} \left(1.0 - \frac{t}{t_{\text{demag}}}\right) \cos(2\pi f_{\text{demag}} t)$$
  $$B_y(t) = B_{\text{max}} \left(1.0 - \frac{t}{t_{\text{demag}}}\right) \sin(2\pi f_{\text{demag}} t)$$
* **Why use it**: Rotational demagnetization sweeps the field vector spherically through the plane, forcing all needles through a complete angular sweep. This erases directional bias and provides a truer randomized state ($S \rightarrow 0$) in 2D interacting grids.

### 2. Simulated Thermal Annealing (`--demag anneal`)
* **Mechanism**: The external field is swept down linearly while a decaying random thermal noise torque (white noise) is applied to all needles.
* **Why use it**: Allows the system configurations to climb out of shallow metastable energy valleys, finding a more isotropic, globally randomized ground state.

---

## 5. System Requirements & Acceleration

`compass_sim` is designed to scale dynamically from small testing runs to large arrays containing thousands of elements.

### CPU Execution (Default)
Uses **NumPy** for vectorization. Best suited for small arrays (e.g. $8 \times 8$ or $12 \times 12$).

### GPU Execution (`--gpu 1`)
Uses **CuPy** to run matrix operations directly on an NVIDIA GPU using custom CUDA kernels.
* **Requirements**: An NVIDIA GPU with CUDA drivers installed, and CuPy (`pip install cupy-cuda12x[ctk]`).
* **Performance Benefit**: Accelerates execution speed significantly for large lattices ($30 \times 30$ and higher) by computing dipole-dipole interactions in parallel.

---

## 6. Command Line Parameters Reference

The `compass_sim` simulator provides a highly configurable interface. The parameters are categorized into six functional groups:

### Lattice Geometry Options
These options configure the spatial dimensions and layout of the needle array:
* `--N <int>`: Number of rows in the lattice grid (default: `8`).
* `--M <int>`: Number of columns in the lattice grid (default: `8`).
* `--R <float>`: Radius of the circumscribed circle enclosing each needle in meters (default: `0.0075`). Center-to-center separation distance between nearest neighbor needles is $2R$.
* `--needle_frac <float>`: Fractional needle length relative to cell diameter $2R$, bounded in $[0.0, 0.8]$. The physical needle length is $L = 2R \cdot \text{needle\_frac}$ (default: `0.6667` = 2/3). Values outside this range are clamped.
* `--geometry {square,triangular,honeycomb}`: The physical layout structure of the grid (default: `square`).

### Physics & Material Properties
These parameters dictate the mechanical and magnetic properties of individual needles, integration time step constraints, and boundary conditions:
* `--moment <float>`: Magnetic moment $m$ of each needle in $\text{A}\cdot\text{m}^2$. If omitted, calculated from geometry assuming saturated steel.
* `--inertia <float>`: Physical moment of inertia $I$ in $\text{kg}\cdot\text{m}^2$. If omitted, calculated assuming a uniform flat rhombic sheet (rhombus) with a central cylindrical pivot hole and including the brass pivot.
* `--pivot_radius <float>`: Radius of the cylindrical pivot at the center of the needle in meters (default: `0.001` = 1.0 mm).
* `--pivot_thickness <float>`: Thickness/height of the cylindrical pivot in meters (default: `0.002` = 2.0 mm).
* `--pivot_density <float>`: Density of the pivot material in $\text{kg/m}^3$ (default: `8500.0` for brass).
* `--pivot_mass <float>`: Explicit mass of the pivot in kg. If provided, overrides density and dimension calculations.
* `--needle_thickness <float>`: Thickness $d$ of the needle sheet in meters (default: `0.00026` = $0.26$ mm).
* `--steel_density <float>`: Density of the needle material in $\text{kg/m}^3$ (default: `7850.0`).
* `--steel_Bsat <float>`: Saturation magnetic flux density $B_{\text{sat}}$ in Tesla (default: `2.0`). Used for calculating magnetic moment ($M_s = B_{\text{sat}} / \mu_0$).
* `--damping <float>`: Viscous damping coefficient $b$ in $\text{N}\cdot\text{m}\cdot\text{s/rad}$ (default: `5.00e-08`).
* `--damping_noise <float>`: Relative variation amplitude for per-needle damping. Damping is randomized as $b_i = \max(0, b \cdot (1 + \text{damping\_noise} \cdot U[-1, 1]))$ (default: `0.0`).
* `--dt_factor <float>`: Fraction of the natural oscillation period $T_0$ to use as the integration time step $\Delta t$ (default: `0.05`).
* `--noise <float>`: Angular noise amplitude in radians applied to initial needle orientations (default: `1.5`).
* `--seed <int>`: Seed for the pseudorandom number generator. If omitted, uses current system nanoseconds.
* `--pbc {0,1}`: Toggle periodic boundary conditions (`0` = disabled, `1` = enabled, default: `0`).
* `--pbc_images <int>`: Number of image shells used to compute long-range periodic dipolar terms (default: `1`).
* `--gpu {0,1}`: Enable GPU acceleration via CuPy (default: `0`).
* `--progress_bar {0,1}`: Toggle the visual terminal progress bar (default: `1`).
* `--halo_mode {order,domains}`: Coloring scheme for needle visualization: `order` shades based on individual deviation, `domains` colors by localized orientation cluster (default: `order`).
* `--domain_tol <float>`: Tolerance angle in degrees used to group needles into discrete domains (default: `15.0`).

### External Magnetic Field Options
These parameters define the field components and profiles:
* `--field_mode {static,hysteresis,sine,pulse,step_pos,step_neg,forc}`: The time-varying protocol for the applied external field $\vec{B}^{\text{ext}}(t)$ (default: `static`).
* `--B_ext <float>`: Magnitude of the external magnetic field in Tesla (default: `0.0`).
* `--phi_ext <float>`: Angle of the external magnetic field in degrees in the $xy$-plane, relative to the $+x$ axis (default: `0.0`).
* `--ext_Bx <float>`: Explicit $x$-component of the external magnetic field in Tesla. Overrides `--B_ext` and `--phi_ext`.
* `--ext_By <float>`: Explicit $y$-component of the external magnetic field in Tesla. Overrides `--B_ext` and `--phi_ext`.
* `--field_delay <float>`: Time delay in seconds before the active field excitation protocol starts (default: `0.0`).
* `--field_freq <float>`: Frequency $f$ in Hz for the sinusoidal external field (default: `1.0`).
* `--t_pulse <float>`: Duration of the field pulse in seconds. Only relevant when `--field_mode` is `pulse`.
* `--torque_tol <float>`: Relative torque threshold used to determine if the lattice has reached mechanical equilibrium (default: `1e-3`).

### FORC & Demagnetization Parameters
Parameters specific to advanced field sweeps and pre-simulation preparation:
* `--forc_n_curves <int>`: Number of minor reversal curves to simulate in FORC mode (default: `30`).
* `--forc_Br_min <float>`: Minimum reversal field $B_{r,\min}$ in Tesla. Defaults to $-B_{\max}$ if unspecified.
* `--forc_t_sat <float>`: Duration in seconds the system is held in positive saturation $+B_{\max}$ prior to starting each curve (default: `0.05`).
* `--forc_t_ramp_down <float>`: Time allocated to ramp from positive saturation to the reversal field $B_r$ (default: `0.10`).
* `--forc_t_ramp_up <float>`: Time allocated to sweep back up from $B_r$ to $+B_{\max}$ (default: `0.20`).
* `--forc_sweep {both,up,down}`: Trace increasing field curves, decreasing field curves, or both (default: `both`).
* `--forc_rate <float>`: Target constant sweep rate $R$ in Tesla/second. If specified, sweep durations are scaled dynamically to enforce a constant rate $dB/dt = R$ across all loops (default: `None`).
* `--demag {off,on,rotational,anneal}`: Selects a pre-relaxation protocol to demagnetize the lattice (default: `off`).
* `--demag_freq <float>`: Rotation/oscillation frequency in Hz of the demagnetization field (default: `2.0`).
* `--demag_cycles <int>`: Number of full decay cycles during demagnetization (default: `20`).
* `--demag_temp <float>`: Noise multiplier factor representing initial thermal energy during simulated annealing (default: `1.0`).
* `--demag_delay <float>`: Settle time in seconds after demagnetization completes (default: `0.0`).
* `--hyst_spacing {linear,log}`: Spacing distribution style of field steps in hysteresis sweeps (default: `linear`).
* `--hyst_log_k <float>`: Concentration parameter controlling spacing density for logarithmic hysteresis sweeps (default: `5.0`).

### Execution & Stop Controls
* `--t_sim <float>`: Total physical time in seconds allocated for the simulation run (default: `2.0`).
* `--t_sim_full {0,1}`: If set to `1`, forces the simulator to run for the full duration of `--t_sim`, bypassing all early physical stopping triggers. Useful for maintaining uniform time-series datasets (default: `0`).

### Visualization & Output Options
* `--video <path>`: Filename or absolute path to export the simulation video (compiled in MP4 format).
* `--frame_every <int>`: Saves one video frame every $N$ integration steps (default: `5`).
* `--make_images {0,1}`: Enable or disable generating summary PNG plots at the end of the simulation (default: `1`).
* `--fps <int>`: Frame rate (frames per second) of the exported MP4 video (default: `24`).
* `--dpi <int>`: Dot-per-inch resolution of saved frames (default: `120`).
* `--keep_frames`: Keep the temporary directory of PNG frames after compilation (default: `False`).
* `--hyst_autoscale {0,1}`: Enable or disable autoscale on the axes of the hysteresis plots (default: `1`).
* `--csv_order {t,B}`: Order of the export columns in the time-series CSV file: `t` writes time first, `B` writes field first (default: `t`).

---

## 7. Output Files & Data Formats

The simulator generates standard files to facilitate analysis and visualization:

### Visual Plots (PNG)
When `--make_images 1` is active, the following summary plots are generated:
* `compass_initial.png`: Shows the initial angular orientation configuration of the needles.
* `compass_equilibrium.png`: Displays the final, relaxed configuration of the needle lattice.
* `compass_comparison.png`: Displays a side-by-side comparison of the initial and final states, along with histograms of angular displacement and local domain size distributions.
* `compass_order_param.png`: Plots the temporal evolution of the lattice-wide order parameter $S(t) = \frac{1}{K} | \sum_{j=1}^K e^{i\theta_j} |$.
* `hysteresis_loop.png` or `forc_curves.png`: Generated during field sweeps. Plots the projected lattice magnetization $M_{\text{proj}}$ as a function of the external field $B_{\text{ext}}$.
* `sine_field.png`: Time-series of both the sinusoidal external field and the lattice's projected magnetization.

### Data Logs (CSV)
The simulator logs numerical data to CSV format. The filename is determined by the field mode and options:
* **General Log (`compass_field_log.csv` or `<video_name>.csv`)**:
  By default, column structure is:
  ```csv
  t_s, B_applied_T, M_proj, S
  ```
  If `--csv_order B` is set, the structure changes to:
  ```csv
  B_applied_T, M_proj, S, t_s
  ```
* **Hysteresis Loop (`hysteresis_loop.csv`)**:
  Logs the sweep trajectory. Columns: `t_s`, `B_T`, `M_proj`, `S`.
* **FORC Sweeps (`forc_curve.csv`)**:
  Logs families of minor loops. Includes the reversal field $B_r$:
  ```csv
  t_s, B_T, M_proj, S, B_r
  ```
* **Sinusoidal Sweeps (`sine_field.csv`)**:
  Logs the sinusoidal response. Columns: `t_s`, `B_T`, `M_proj`, `S`.

### Simulation Videos (MP4)
If `--video <name>` is specified, a high-quality video is rendered showing the dynamic spatial evolution of the lattice. Needle colors dynamically correspond to their instantaneous angle or local domain affinity depending on the coloring mode.

---

## 8. Step-by-Step Practical Examples

Here are common recipes for simulating different states.

### Example 1: Basic Equilibrium Relax (No Field)
Run the default $8 \times 8$ square lattice from a highly randomized state and let it settle under dipolar interaction into an ordered ground state.
```bash
python3 compass.py --t_sim 2.0 --noise 2.5
```
**Output**: Generates:
* `compass_initial.png`: Highly randomized needle orientations.
* `compass_equilibrium.png`: The final state (typically micro-domains or antiferromagnetic configurations).
* `compass_order_param.png`: Evolution of the order parameter $S(t)$ over time.

### Example 2: Hysteresis Loop
Trace a complete magnetization cycle ($M$ vs. $B_{\text{ext}}$) for a $10 \times 10$ triangular lattice with a maximum applied field of $1.5\text{ mT}$ at $30^{\circ}$ offset.
```bash
python3 compass.py --geometry triangular --N 10 --M 10 \
    --field_mode hysteresis --B_ext 1.5e-3 --phi_ext 30 \
    --t_sim 5.0 --damping 5e-7
```
**Output**: Displays the hysteresis loop plot `compass_comparison.png` and exports data to `hysteresis_loop.csv`.

### Example 3: Simulating Demagnetization
Run a pulse simulation where the lattice is first rotational-demagnetized at $5\text{ Hz}$ over $40$ cycles, allowed to rest for $0.2\text{ s}$, then subjected to a $100\ \mu\text{T}$ pulse.
```bash
python3 compass.py --field_mode pulse --B_ext 100e-6 --phi_ext 90 \
    --demag rotational --demag_freq 5.0 --demag_cycles 40 --demag_delay 0.2 \
    --t_pulse 0.5 --t_sim 1.5 --video pulse_demag.mp4
```
**Output**: Saves an MP4 video `pulse_demag.mp4` rendering the initial biaxial demagnetization sweep, the rest period, and the subsequent pulse response.

### Example 4: Large-Scale Run on GPU
Run a large $30 \times 30$ square lattice simulation (900 needles) using CUDA acceleration.
```bash
python3 compass.py --N 30 --M 30 --gpu 1 --t_sim 2.5 --video large_grid.mp4
```

---

---

## 9. Running Sweep Campaigns

For parametric analysis (e.g. studying how damping parameters alter hysteresis loops), use the campaign wrapper script:

```bash
python3 sweep_damping_hysteresis.py --out_dir ./damping_campaign
```

### Script Execution Parameters:
* `--out_dir <dir>`: Target directory for CSV field logs and JSON metadata.
* `--resume`: Scans the target folder and resumes progress, skipping already completed parameter combinations.
* `--dry_run`: Outputs the list of scheduled simulation parameters without running them.
* `--gpu 1`: Enables GPU acceleration inside the sweep campaign runs.

---

## 10. Summary Reference Table of Parameters

| Parameter | Type / Options | Default | Brief Description |
| :--- | :--- | :--- | :--- |
| **Lattice Geometry Options** | | | |
| `--N` | `int` | `8` | Number of rows in the lattice grid. |
| `--M` | `int` | `8` | Number of columns in the lattice grid. |
| `--R` | `float` | `0.0075` | Cell circumscribed radius (m); center spacing is $2R$. |
| `--needle_frac` | `float` | `0.6667` | Needle size ratio relative to $2R$, range $[0.0, 0.8]$. |
| `--geometry` | `square`, `triangular`, `honeycomb` | `square` | Physical spatial layout arrangement. |
| **Physics & Material Properties** | | | |
| `--moment` | `float` | `None` | Needle magnetic moment $m$ ($A\cdot m^2$). Calculated if omitted. |
| `--inertia` | `float` | `None` | Needle moment of inertia $I$ ($kg\cdot m^2$). Calculated if omitted. |
| `--pivot_radius` | `float` | `0.001` | Radius of the cylindrical pivot (m). |
| `--pivot_thickness` | `float` | `0.002` | Thickness/height of the central cylindrical pivot (m). |
| `--pivot_density` | `float` | `8500.0` | Density of pivot material ($kg/m^3$). |
| `--pivot_mass` | `float` | `None` | Explicit mass of the central pivot (kg). Overrides density/dimensions. |
| `--needle_thickness` | `float` | `0.00026` | Thickness $d$ (m) used for moment/inertia auto-calibration. |
| `--steel_density` | `float` | `7850.0` | Density of needle material ($kg/m^3$). |
| `--steel_Bsat` | `float` | `None` | Saturation flux density $B_{\text{sat}}$ (T). Default uses $2.0$ T. |
| `--damping` | `float` | `5.00e-08` | Viscous damping coefficient $b$ ($N\cdot m\cdot s/rad$). |
| `--damping_noise` | `float` | `0.0` | Relative variation factor in damping coefficient per needle. |
| `--dt_factor` | `float` | `0.05` | Integration step factor $\Delta t / T_0$ for Velocity-Verlet. |
| `--noise` | `float` | `1.5` | Angular noise amplitude (rad) for initial orientations. |
| `--seed` | `int` | `None` | Pseudorandom number generator seed. Uses time-ns if omitted. |
| `--pbc` | `0`, `1` | `0` | Toggle periodic boundary conditions. |
| `--pbc_images` | `int` | `1` | Periodic replication shell summation depth. |
| `--gpu` | `0`, `1` | `0` | Toggle GPU acceleration via CuPy. |
| `--progress_bar` | `0`, `1` | `1` | Toggle visual interactive CLI progress bar. |
| `--halo_mode` | `order`, `domains` | `order` | Halo coloring mode for needle rendering. |
| `--domain_tol` | `float` | `15.0` | Angle threshold (degrees) for cluster-grouping domains. |
| **External Magnetic Field Options** | | | |
| `--field_mode` | `static`, `hysteresis`, `sine`, `pulse`, `step_pos`, `step_neg`, `forc` | `static` | The excitation profile for the applied field $\vec{B}^{\text{ext}}(t)$. |
| `--B_ext` | `float` | `0.0` | External field magnitude limit (T). |
| `--phi_ext` | `float` | `0.0` | Applied field angle (degrees) relative to $+x$. |
| `--ext_Bx` | `float` | `None` | Direct $x$-component of the external magnetic field (T). |
| `--ext_By` | `float` | `None` | Direct $y$-component of the external magnetic field (T). |
| `--field_delay` | `float` | `0.0` | Delay duration (s) before applying field protocol. |
| `--field_freq` | `float` | `1.0` | Sinusoidal AC field frequency (Hz). |
| `--t_pulse` | `float` | `None` | Pulse duration (s) under pulse mode. |
| `--torque_tol` | `float` | `1e-3` | Rotational equilibrium early stop tolerance. |
| **FORC & Demagnetization Parameters** | | | |
| `--forc_n_curves` | `int` | `30` | Number of First-Order Reversal Curves to trace. |
| `--forc_Br_min` | `float` | `None` | Minimum reversal field limit (T). Defaults to $-B_{\max}$. |
| `--forc_t_sat` | `float` | `0.05` | Holding time at positive saturation field (s). |
| `--forc_t_ramp_down` | `float` | `0.10` | Holding/sweep ramp down time to $B_r$ (s). |
| `--forc_t_ramp_up` | `float` | `0.20` | Holding/sweep ramp up time back to saturation (s). |
| `--forc_sweep` | `both`, `up`, `down` | `both` | Curve sweep directions to save and record in CSV. |
| `--forc_rate` | `float` | `None` | Constant $dB/dt$ rate (T/s) for dynamic curve duration scaling. |
| `--demag` | `off`, `on`, `rotational`, `anneal` | `off` | Selected demagnetization method prior to simulation. |
| `--demag_freq` | `float` | `2.0` | Frequency (Hz) of demagnetization fields. |
| `--demag_cycles` | `int` | `20` | Number of field decay cycles for demagnetization. |
| `--demag_temp` | `float` | `1.0` | Stochastic noise factor for simulated annealing. |
| `--demag_delay` | `float` | `0.0` | Rest settle time (s) after demagnetization ends. |
| `--hyst_spacing` | `linear`, `log` | `linear` | Spacing distribution type for hysteresis field steps. |
| `--hyst_log_k` | `float` | `5.0` | Log-spacing density parameter for hysteresis. |
| **Execution & Stop Controls** | | | |
| `--t_sim` | `float` | `2.0` | Total simulation physical time (s). |
| `--t_sim_full` | `0`, `1` | `0` | Enforce full $t_{\text{sim}}$ duration execution without early stops. |
| **Visualization & Output Options** | | | |
| `--video` | `str` | `None` | Filepath to save rendered MP4 animation. |
| `--frame_every` | `int` | `5` | Save frame image every $N$ integration steps. |
| `--make_images` | `0`, `1` | `1` | Enable summary and comparison PNG plot exports. |
| `--fps` | `int` | `24` | Frames per second of compiled video. |
| `--dpi` | `int` | `120` | Dot-per-inch resolution of frames. |
| `--keep_frames` | flag | `False` | Keep raw PNG frames folder after video compile. |
| `--hyst_autoscale` | `0`, `1` | `1` | Enable autoscale on the axes of hysteresis plots. |
| `--csv_order` | `t`, `B` | `t` | Column structure in field log CSV file (t-first or B-first). |
