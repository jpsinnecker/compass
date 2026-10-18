# Walkthrough - Default Physical Rhombic Needle Geometry (V79)

I have successfully updated the simulator and all accompanying scripts and documentation to use the new default physical parameters derived from real physical needle dimensions (10 mm long, 3 mm wide, 0.26 mm thick, with a 2.0 mm diameter x 2.0 mm height brass pivot pin).

## Physical Derivations

For a rhombic needle of circumscribed radius $R=7.5$~mm with `needle_frac = 2/3`:
- **Length $L$**: $2R \cdot \text{needle\_frac} = 15 \cdot (2/3) = 10$~mm.
- **Width $W$**: $10 \cdot 0.30 = 3$~mm.
- **Hole Radius $r_p$**: $1.0$~mm.
- **Hole Thickness/Height $h_p$**: $2.0$~mm.
- **Net Blade Volume**: $V_{\text{blade}} = (\text{Area}_{\text{rhombus}} - \text{Area}_{\text{hole}}) \cdot d \approx (1.5\times10^{-5} - \pi\times10^{-6}) \cdot 0.26\times10^{-3} \approx 3.083\times10^{-9}$~m³.
- **Blade Mass**: $M_{\text{blade}} = \rho_{\text{steel}} \cdot V_{\text{blade}} \approx 7850 \cdot 3.083\times10^{-9} \approx 2.42\times10^{-5}$~kg = $0.0242$~g.
- **Pivot Mass**: $M_{\text{pivot}} = \rho_{\text{brass}} \cdot (\pi\cdot r_p^2 \cdot h_p) \approx 8500 \cdot \pi\cdot10^{-6}\cdot0.002 \approx 5.34\times10^{-5}$~kg = $0.0534$~g.
- **Total Mass**: $M_{\text{total}} \approx 7.76\times10^{-5}$~kg = $0.0776$~g.
- **Blade Inertia $I_{\text{blade}}$**: Subtracting the pivot hole from the solid rhombus:
  $$I_{\text{solid}} = \frac{1}{24} M_{\text{solid}} (L^2 + W^2) \approx \frac{1}{24} (3.06\times10^{-5}) (0.01^2 + 0.003^2) \approx 1.39\times10^{-10}\text{ kg}\cdot\text{m}^2$$
  $$I_{\text{hole}} = \frac{1}{2} M_{\text{hole}} r_p^2 \approx \frac{1}{2} (6.41\times10^{-6}) (0.001^2) \approx 3.2\times10^{-12}\text{ kg}\cdot\text{m}^2$$
  $$I_{\text{blade}} = I_{\text{solid}} - I_{\text{hole}} \approx 1.358\times10^{-10}\text{ kg}\cdot\text{m}^2$$
- **Pivot Inertia $I_{\text{pivot}}$**: $0.5 \cdot M_{\text{pivot}} \cdot r_p^2 \approx 0.5 \cdot (5.34\times10^{-5}) (0.001^2) \approx 2.67\times10^{-11}\text{ kg}\cdot\text{m}^2$.
- **Total Inertia**: $I_{\text{total}} = I_{\text{blade}} + I_{\text{pivot}} \approx 1.625\times10^{-10}\text{ kg}\cdot\text{m}^2$.
- **Magnetic Moment $m$**: $M_s \cdot V_{\text{blade}} \approx 1.59\times10^6 \cdot 3.083\times10^{-9} \approx 4.90\times10^{-3}\text{ A}\cdot\text{m}^2$.

## Changes Made

### 1. Main Simulator Logic
- **[compass.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/compass.py)**
  - Updated default values: circumscribed radius $R = 0.0075$ ($7.5$ mm), `needle_frac = 2.0/3.0` ($66.7\%$), steel thickness $d = 0.00026$ ($0.26$ mm), pivot radius $r_p = 1.0\times10^{-3}$ ($1.0$ mm), pivot thickness $h_p = 2.0\times10^{-3}$ ($2.0$ mm).
  - Updated width-to-length ratio from $0.22$ to $0.30$.
  - Subtracted the central pivot hole area/volume from the blade mass, moment of inertia, and magnetic moment auto-calibration.
  - Corrected blade mass formula in the banner printout to correctly show the rhombus blade mass (minus the hole).
  - Bumped code version to `V79` in the banner and comments.

### 2. Campaign & Analysis Scripts
- **[damping_sweep.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/damping_sweep.py)** & **[damping_sweepV02.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/damping_sweepV02.py)**
  - Updated argument parser defaults (`--R 0.0075`, `--needle_frac 2.0/3.0`, `--needle_thickness 0.26e-3`, `--pivot_radius 1.0e-3`, `--pivot_thickness 2.0e-3`).
  - Updated default constants (`R_DEFAULT = 0.0075`, `NEEDLE_FRAC = 2.0/3.0`).
  - Updated `derive_needle_geometry` to use width ratio $0.30$ and thickness $0.26$ mm.
  - Forwarded `pivot_radius` to `compute_moment_from_geometry` inside `compute_physical_params`.
- **[sweep_damping_hysteresis.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/sweep_damping_hysteresis.py)**
  - Updated default constants (`R_DEFAULT = 0.0075`, `NEEDLE_FRAC = 2.0/3.0`, `NEEDLE_THICKNESS = 0.26e-3`).
  - Updated `_compute_derived_params` to use width ratio $0.30$ and default pivot dimensions.

### 3. Project Documentation
- **[USER_GUIDE.md](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/USER_GUIDE.md)**
  - Updated parameter descriptions and default values in the documentation text and summary reference table.
- **[project_report.tex](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/project_report.tex)**
  - Updated defaults in the parameters table (Table 1) and introductory physics sections.
  - Bumped report version to `V50` and added entries for code `V79` and report `V50` to the history tables.
  - Compiled the updated LaTeX report successfully.

## Verification & Test Results

We ran a verification simulation in the scratch directory:
```bash
python3 compass.py --N 4 --M 4 --t_sim 0.01
```
- Console output correctly prints the new version, geometry, and computed physical variables:
  ```text
  ══════════════════════════════════════════════════════════════
    Version      : V79
    Developer    : J.P. Sinnecker
    Backend      : CPU (NumPy)  [--gpu 0]
    Array        : square  4x4 needles
    PBC          : off
    Radius R     : 0.75 cm  (2R = 1.50 cm)
    Needle R     : 67% of 2R = 1.00 cm
    Mag Moment.  : 0.004902 A.m2  per needle
    Inertia      : 1.625e-10 kg.m2  per needle (blade: 2.42e-02g, pivot: 5.34e-02g)
    Damping      : 5.000e-08 N.m.s/rad  (air)
    Q            : 3.043e-01  (overdamped)
  ══════════════════════════════════════════════════════════════
  ```
- The simulation runs and logs correctly with the new physical parameter defaults.
- The compiled `project_report.pdf` has been fully generated.
