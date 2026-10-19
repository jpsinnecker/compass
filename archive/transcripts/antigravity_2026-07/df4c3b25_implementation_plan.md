# Implementation Plan — Bidirectional Field Box Bar (V66) and Report Update (V33)

This plan outlines the changes required to update the external field intensity bar in `compass.py` to grow bidirectionally (left for negative fields, right for positive fields, and zero at zero field, matching the provided design). It also details renaming the project report to `project_report.tex`, incrementing the report version to `V33`, and documenting the version changes.

## Proposed Changes

### 1. Python Code

#### [MODIFY] [compass.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/compass.py)

- **Change code version to V66**:
  - Update line 3 comment to reflect version `V66`.
  - Update line 2208 `_print(f"  Version      : V66")`.
- **Implement bidirectional field bar**:
  - Locate `draw_ext_field_on_lattice` function.
  - Replace the single-direction progress bar rendering with a centered bidirectional rendering:
    - Compute the signed fraction `f = b_val / B_max_ref`.
    - If `f >= 0`, draw the active rectangle starting at `center_x = px + bar_w / 2.0` extending rightwards by `(bar_w / 2.0) * f`.
    - If `f < 0`, draw the active rectangle starting at `center_x - (bar_w / 2.0) * abs(f)` extending rightwards by `(bar_w / 2.0) * abs(f)`.
    - Draw a solid red vertical reference line in the center at `center_x`.
    - Retain the background rectangle and outer borders as they are.

#### [NEW] [freeze_v66.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/freeze_v66.py)
- Create `freeze_v66.py` to copy `compass.py` to `compass_simV66.py`.

---

### 2. Project Report

#### [NEW] [project_report.tex](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/project_report.tex)
- Copy the content of `relatorio_projeto_V22_EN.tex` to `project_report.tex`.
- In `project_report.tex`:
  - Update preamble titles, version notes, and compilation comments to refer to `project_report.tex` and `Report version V33`.
  - Update abstract and text to cover history up to code version `V66` (`code-version history (V20--V66)`).
  - Add a new row to the version log table:
    ```latex
    V33 & 06/07/2026 & Documented code V66 (HUD external field intensity bar redesigned to grow bidirectionally from center: left for negative fields, right for positive fields, zero at zero field, with a central red zero-reference line). \\
    ```

#### [MODIFY] [build_report.sh](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/build_report.sh)
- Update `TEX` parameter to `"project_report"` so that running the build script automatically compiles `project_report.tex`.

---

## Verification Plan

### Automated Tests
- Run Python syntax checks on `compass.py` and test run to verify compilation correctness:
  ```bash
  python3 -m py_compile compass.py
  ```
- Run a short check by running `compass.py` to generate the default plots:
  ```bash
  python3 compass.py --B_ext 50e-6 --phi_ext 45 --t_sim 0.01 --make_images 1
  ```
- Compile the report to PDF:
  ```bash
  bash build_report.sh
  ```

### Manual Verification
- Inspect the generated plots/figures to ensure that:
  - The external field box shows the yellow bar growing from the center to the left/right depending on the field sign.
  - A red line is shown in the center of the bar.
  - At zero field, the bar is empty and the red line is at the center.
- Verify the generated PDF report `project_report.pdf` is correctly updated.
