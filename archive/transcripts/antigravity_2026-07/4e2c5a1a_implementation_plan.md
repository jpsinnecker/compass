# FORC Analysis — Integration into the Damping Sweep Study

## Background

`compass.py` already implements `field_mode='forc'` in its `relax()` function.
The FORC engine:
- Runs **n** reversal curves, each: saturate → ramp down to B_r → ramp back to B_max
- Logs rows of `[t, B, M, S, B_r, sweep_dir, curve_idx]`
- Key parameters: `forc_n_curves`, `forc_Br_min`, `forc_t_sat`, `forc_t_ramp_down`, `forc_t_ramp_up`, `forc_sweep`, `forc_rate`

The FORC distribution ρ(B_r, B) = -½ ∂²M/∂B_r∂B is computed in **analysis**, not inside the sweep.

## Open Questions

> [!IMPORTANT]
> **Q1 — Scope of FORC sweep**: Should FORC run for every (geometry, damping, seed) triple (same as hysteresis), or only for a **selected subset** of Q values (e.g., under-damped, critical, over-damped) to keep runtime reasonable? FORC with 30 curves costs ~6× more than a single hysteresis loop per run.

> [!IMPORTANT]
> **Q2 — forc_rate vs. fixed times**: Should ramp times scale proportionally with sweep depth (`forc_rate` = constant T/s), or use fixed `forc_t_ramp_down / forc_t_ramp_up`? Constant rate is more physically consistent; fixed times are simpler.

> [!NOTE]
> **Q3 — B_r min**: The default in compass.py is `forc_Br_min = -B_max`. Is this the right range, or should it stop at 0 (only the irreversible region)?

## Proposed Changes

---

### `damping_sweepV02.py`

#### [MODIFY] [damping_sweepV02.py](file:///Users/jps/Dropbox/simula%20agulhas/damping_sweepV02.py)

1. **New function `run_one_forc()`** — mirrors `run_one_hysteresis()` but calls
   `cs.relax(..., field_mode='forc', forc_n_curves=..., forc_Br_min=..., forc_rate=..., ...)`.
   Returns `(field_log, meta)` where `field_log` columns are `[t, B, M, S, B_r, sweep_dir, curve_idx]`.

2. **New `save_forc_csv()`** — saves the 7-column FORC log (vs. the 4-column hysteresis log).

3. **`run_job()` extension** — adds `stage == 'forc'` branch.

4. **`main()` additions**:
   - New `forc_dir = out_dir/forc/` output directory.
   - New CLI flags:
     - `--skip_forc_stage 0` (default: run FORC)
     - `--forc_n_curves 20` (default for production; 3 for quick_test)
     - `--forc_rate` (T/s; None → fixed times)
     - `--forc_Br_min` (default: -B_max)
     - `--forc_t_sat / forc_t_ramp_down / forc_t_ramp_up` (all defaulting to compass.py defaults)
   - Stage 3 job list added to `jobs[]`.
   - `quick_test` override: `forc_n_curves = 3`.
   - `MANIFEST_FIELDS` gains `forc_n_curves` for FORC rows.
   - `common` dict gains `forc_dir` and all forc_ parameters.

---

### `damping_sweep.py` (serial V01)

#### [MODIFY] [damping_sweep.py](file:///Users/jps/Dropbox/simula%20agulhas/damping_sweep.py)

Same additions as V02 but in the simpler serial loop structure (no `run_job()` pattern).

---

### `damping_sweep_analysis.py`

#### [MODIFY] [damping_sweep_analysis.py](file:///Users/jps/Dropbox/simula%20agulhas/damping_sweep_analysis.py)

1. **`load_forc_data(row)`** — reads a FORC CSV, returns a DataFrame with columns
   `[t, B, M, S, B_r, sweep_dir, curve_idx]`.

2. **`compute_forc_distribution(forc_df, n_grid=150, smooth_sigma=2)`** — computes the FORC
   distribution ρ(B_r, B) via:
   - For each (B_r, B) point, interpolate M onto a regular 2D grid
   - Apply Gaussian smoothing (σ = `smooth_sigma` grid cells)
   - Compute mixed second derivative: ρ = -½ ∂²M/∂B_r∂B  (finite differences)

3. **`plot_forc_diagram(rho, Br_grid, B_grid, out_path, title)`** — contour plot with:
   - Horizontal axis: B (applied field during ascending branch) in mT
   - Vertical axis: B_r (reversal field) in mT
   - Colour scale: blue → white → red, symmetric around 0
   - Coercivity contour (B_c = (B - B_r)/2 axis) overlaid as dashed line
   - Interaction field contour (B_u = (B + B_r)/2 axis)

4. **`process_forc_run(row)`** — calls 1 + 2 for one manifest row; returns `(rho, Br_grid, B_grid, meta)`.

5. **`plot_forc_vs_Q()`** — panel plot: one FORC diagram per Q value, one panel per geometry.
   Output: `forc_diagrams.png`

6. **`main()` additions**:
   - Read `forc` rows from manifest.
   - New Stage 3 loop (with progress bar `label="FORC"`).
   - Aggregate FORC distributions per geometry across seeds → mean ρ per (geometry, Q).
   - Call `plot_forc_vs_Q()`.
   - Add `forc_diagrams.png` to the "Files written" list.
   - Remove the `--forc_dir` stub / `Ignored` print.

---

## Verification Plan

### Automated Tests
```bash
# Quick pipeline end-to-end (3 curves, tiny grid):
python3 damping_sweepV02.py --out_dir /tmp/forc_test --quick_test 1 \
    --n_workers 2 --geometries square --forc_n_curves 3

python3 damping_sweep_analysis.py \
    --sweep_dir /tmp/forc_test --out_dir /tmp/forc_test/analysis
```

### Manual Verification
- Confirm `forc_diagrams.png` is generated without error.
- Check that the FORC distribution has the expected elongated ridge along the B_c axis for overdamped Q and a broader distribution for underdamped Q.
- Confirm the FORC CSV has 7 columns and B_r values span the expected range.
