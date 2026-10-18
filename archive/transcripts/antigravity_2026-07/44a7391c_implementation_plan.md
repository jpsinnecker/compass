# Implementation Plan - Default Physical Rhombic Needle Geometry (V79)

This plan details the changes required to set the default needle geometry parameters in `compass.py` to match the real physical needles: 10 mm long, 3 mm wide, 0.26 mm steel thickness, with a central brass cylinder pivot (2.0 mm diameter x 2.0 mm height).

## Proposed Physical Model & Defaults

1. **Rhombic Needle Geometry**:
   - The needle blade is modeled as a **rhombus (diamond)** rather than a rectangle.
   - **Solid Rhombus Area**: $A_{\text{solid}} = \frac{1}{2} L W$
   - **Solid Rhombus Inertia**: $I_{\text{solid}} = \frac{1}{24} M_{\text{solid}} (L^2 + W^2)$
   
2. **Central Pivot & Hole Subtraction**:
   - A central cylindrical hole of radius $R_p = 1.0$ mm is subtracted from the steel blade's mass, volume, and moment of inertia.
   - A central brass cylinder pivot of radius $R_p = 1.0$ mm and height $h_p = 2.0$ mm (density $8500.0$ kg/m³) is added to the needle's mass and moment of inertia.

3. **New Defaults in `compass.py`**:
   - Fallback `MOMENT_DEFAULT = 4.902e-3` A·m²
   - Fallback `INERTIA_DEFAULT = 1.625e-10` kg·m²
   - `--R` default: `0.0075` m (7.5 mm)
   - `--needle_frac` default: `2.0/3.0` (0.6667)
   - `--needle_thickness` default: `0.26e-3` m (0.26 mm)
   - `--pivot_radius` default: `1.0e-3` m (1.0 mm)
   - `--pivot_thickness` default: `2.0e-3` m (2.0 mm)
   - Width ratio in `main()`: `needle_width = needle_len * 0.30`

---

## Proposed Changes

### 1. Simulator Core

#### [MODIFY] [compass.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/compass.py)
- Update `MOMENT_DEFAULT`, `INERTIA_DEFAULT`, and `NEEDLE_THICKNESS_DEFAULT`.
- Modify `compute_inertia_from_geometry` to calculate rhombic blade inertia, subtract the cylindrical hole, and add the brass pivot cylinder.
- Modify `compute_moment_from_geometry` to subtract the pivot hole volume from the magnetic steel blade volume.
- Update argument defaults in the CLI parser for `--R`, `--needle_frac`, `--needle_thickness`, `--pivot_radius`, and `--pivot_thickness`.
- Update `needle_width` calculation in `main()` to use a `0.30` ratio.
- Bump version to `V79` in banner, comments, and reports.

### 2. Satellite Scripts

#### [MODIFY] [damping_sweep.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/damping_sweep.py)
#### [MODIFY] [damping_sweepV02.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/damping_sweepV02.py)
#### [MODIFY] [sweep_damping_hysteresis.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/sweep_damping_hysteresis.py)
- Update default geometry parameters (`R`, `thickness`, `needle_frac`, width ratio `0.30`, and pivot defaults) to match the new defaults of `compass.py` for physical consistency.

### 3. Documentation

#### [MODIFY] [project_report.tex](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/project_report.tex)
- Update defaults table, document the rhombic physical model equations, and log report version `V49` updates.
- Compile updated PDF using `pdflatex`.

#### [MODIFY] [USER_GUIDE.md](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/USER_GUIDE.md)
- Update geometry default listings and descriptions.

---

## Verification Plan

### Automated Tests
- Run inline check command to print the default parameter block:
  ```bash
  python3 compass.py --N 4 --M 4 --t_sim 0.01
  ```
- Verify that default inertia, moment, and dimensions match:
  - `Inertia: 1.625e-10 kg.m2`
  - `Mag Moment: 0.004902 A.m2`
  - `Radius R: 0.75 cm (2R = 1.50 cm)`
  - `Needle R: 67% of 2R = 1.00 cm`
  - `Q: 3.639e-03 (overdamped)` for default damping `5e-8`.
