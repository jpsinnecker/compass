# Update Animation Frames for Hysteresis Mode

This plan modifies the `compass.py` script so that when it is run in `hysteresis` mode, the generated video (or frames) will contain two subplots side-by-side:
- **Left:** The compass needles (the current existing visual).
- **Right:** The hysteresis curve ($M \times B$) being dynamically updated, with a dot marking the current state.

## Proposed Changes

### `compass.py`

#### [MODIFY] compass.py

1. **Update `plot_state` function signature and internals:**
   - Add `ax=None` to the parameters list.
   - If `ax is None`, compute `figsize_inches` and create `fig, ax = plt.subplots(...)` as it currently does.
   - If `ax` is provided, do not create a new figure, simply set `fig = ax.figure`.

2. **Update `_save_frame` function inside `simulate_network`:**
   - Detect if `field_mode == 'hysteresis'`.
   - If so, create a figure with two side-by-side axes: `fig, (ax_net, ax_curve) = plt.subplots(1, 2, figsize=(14, 6), facecolor='#1A1A2E')`.
   - Call `plot_state` passing `ax=ax_net`.
   - Use `ax_curve` to plot the `hysteresis_log` up to the current step.
     - Plot `B_arr` vs `M_arr`.
     - Draw a golden dot at the last point: `(B_arr[-1], M_arr[-1])`.
     - Format the axes to match the aesthetics of the existing `_plot_hysteresis` function (dark background, grid, labels).
     - Fix the limits of the x-axis to $[-B_{max} \times 1.1, B_{max} \times 1.1]$ and y-axis to $[-1.1, 1.1]$ so the plot frame stays stable while the curve is drawn.
   - If not in hysteresis mode, call `plot_state` as usual.

## User Review Required

> [!NOTE]
> The new video format will have a `14x6` aspect ratio (wide) specifically during `field_mode='hysteresis'` to fit the curve nicely on the right. Is this dimension suitable for your needs?

## Verification Plan

### Manual Verification
- Run the simulation in hysteresis mode with a short duration or low number of points and `--video test_hysteresis.mp4`.
- Check if the MP4 video is rendered with the needles on the left and a live-drawing curve on the right.
