# Implementation Plan: Multi-Stage Progress Bar

This plan proposes dividing the simulator's command-line progress bar into two distinct phases when a demagnetization pre-treatment protocol is enabled.

## Proposed Changes

### Terminal Output & Visuals

#### [MODIFY] [compass.py](file:///Users/jps/Dropbox/simula%20agulhas/compass.py)

1. **`_print_progress_bar` parameter extension**:
   - Add a `label="Integrating"` parameter to `_print_progress_bar`.
   - Update `line1` of the progress bar to show `{label}` instead of hardcoded `"Integrating"`.
2. **Phase Partitioning in `relax(...)` loop**:
   - Inside the integration loop, identify if the current step is within the demagnetization phase (`step <= n_steps_demag`) or the subsequent simulation/measurement phase (`step > n_steps_demag`).
   - If in the **Demagnetization phase**:
     - Progress fraction: `step / n_steps_demag`
     - Step count: `step` out of `n_steps_demag`
     - Time: `t_now`
     - Label: `"Demagnetizing"`
   - If in the **Simulation phase**:
     - Progress fraction: `(step - n_steps_demag) / (n_steps - n_steps_demag)`
     - Step count: `step - n_steps_demag` out of `n_steps - n_steps_demag`
     - Time: `t_phys`
     - Label: `"Hysteresis"` (or `"FORC sweep"`, `"Relaxing"`, `"Integrating"` depending on `field_mode`).
3. **Transition Finalization**:
   - At the exact step where demagnetization completes (`step == n_steps_demag`):
     - Force-print the final 100% demagnetization progress bar.
     - Finalize the bar via `_print_progress_bar_finish()`.
     - Print a transition header: `  Demagnetization complete. Starting [simulation phase]...`
     - Reset performance print counters so that the subsequent phase progress bar starts counting cleanly.

---

## Verification Plan

### Automated Tests
Run a simulation with demagnetization enabled to observe the two separate progress bars:
```bash
python3 compass.py --field_mode static --B_ext 0.005 --demag rotational --demag_cycles 5
```

### Manual Verification
- Verify that the terminal displays two distinct progress bars sequentially (first "Demagnetizing" up to 100%, then the transition header, and finally "Relaxing" up to 100%).
