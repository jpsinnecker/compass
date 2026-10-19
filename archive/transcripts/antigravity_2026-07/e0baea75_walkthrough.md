# Walkthrough: Hysteresis Video Layout

The modifications to `compass.py` to support side-by-side video rendering in `hysteresis` mode have been completed. 

## Changes Made

### 1. Extracted Axis Support in `plot_state`
- `plot_state` now accepts an optional `ax` parameter.
- If `ax` is provided, it draws the compass network directly onto that axis instead of creating a standalone figure. This avoids overlapping figures and allows `plot_state` to be embedded in subplots.

### 2. Side-by-Side Layout in `_save_frame`
- Inside the frame generation loop `_save_frame`, the script now checks if `field_mode == 'hysteresis'`.
- If true, a `14x6` figure is created with two subplots:
  - **Left subplot:** The compass network is rendered here using the modified `plot_state(..., ax=ax_net)`.
  - **Right subplot:** The hysteresis curve ($M \times B$) is plotted dynamically using the `hysteresis_log` data up to the current simulation step.
- A golden dot (`o`) is plotted on the curve to indicate the current state of the system in the cycle.
- The axes limits for the hysteresis curve are fixed dynamically based on `B_ext_max` to prevent the view from jumping around as the curve is drawn.
- When `field_mode != 'hysteresis'`, it gracefully defaults back to the previous behavior (a single square frame).

## Validation

> [!WARNING]
> Due to an environment configuration issue on your system (Anaconda's Python is throwing a `ModuleNotFoundError: No module named 'encodings'` when executing python), I was not able to perform an end-to-end rendering test on the terminal.
>
> You can try running the following command in a fresh terminal once your Python environment is fixed to see the result:
> ```bash
> python compass.py --rows 10 --cols 10 --field_mode hysteresis --B_ext 5 --t_sim 0.5 --fps 10 --video hysteresis_demo.mp4
> ```
