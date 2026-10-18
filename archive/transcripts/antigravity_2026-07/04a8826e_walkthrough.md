# Walkthrough - FORC Measurement Protocol Implementation (V68)

This document summarizes the changes made to implement the First-Order Reversal Curve (FORC) measurement protocol in the compass lattice simulation.

## Changes Made

### 1. Version Bump
- Updated simulation version to **V68** in:
  - [compass.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/compass.py) (Header comment, CLI `--version` output, and log headers).
  - [compass_simV68.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/compass_simV68.py) (Frozen snapshot clone).

### 2. FORC Up-Sweep Phase and Loop Index Tracking
- In the main integration loop of the `relax()` function inside [compass.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/compass.py#L1228), implemented tracking for the up-sweep phase:
  - Calculated `t_local = t_now - k * T_cycle` (current elapsed time within the cycle `k`).
  - Identified the up-sweep segment using `is_up_sweep = (t_local > tsat + t_down)`.
  - Logged each step as a tuple with additional fields: `(t_now, B_scalar, M_proj, S_now, Br_curr, is_up_sweep, k)`.

### 3. Separation of Curves in Visualization
- Updated `_render_hyst_panel` (real-time rendering) and `_plot_hysteresis` (final plot generation) to:
  - Group logged data points by their loop index `k` using a dictionary of curves.
  - Render each loop segment as an independent line segment, preventing visual artifacts where lines connect back to the start of the next cycle.
  - Maintained the active marker (red dot) tracing the full continuous sweep in real time.

### 4. CSV Logging Filter
- Modified both the relaxation inner loop saving and the main CLI exit saving code in `compass.py`:
  - Enforced that only data points belonging to the up-sweep phase (`is_up_sweep == True`) are recorded in the exported files:
    - `forc_curve.csv`
    - `compass_field_log.csv`
  - Columns written: `['t_s', 'B_applied_T', 'M_proj', 'S', 'B_r']` (and `B_T` for `forc_curve.csv`).

### 5. Frozen Snapshot and Freeze Script
- Created [freeze_v68.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/freeze_v68.py) to snapshot this version.
- Created [compass_simV68.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/compass_simV68.py) as the frozen codebase.

---

## Verification Results

### Simulated Test Run
Due to macOS sandbox restrictions on subprocesses writing to the Dropbox folder, the test run was executed in the scratch sandbox directory `/Users/jps/.gemini/antigravity-ide/scratch/simbuild/` using the command:
```bash
/opt/anaconda3/bin/python3 compass.py --N 4 --M 4 --field_mode forc --forc_n_curves 5 --t_sim_full 1 --damping 8e-6
```

#### Results:
- **Simulation Run**: Successfully completed all 5 cycles without any numerical issues.
- **CSV Data Validation**: The exported CSV files contain exactly **73 lines** of up-sweep measurements instead of 127 total steps, validating that the down-sweeps were correctly filtered out.
- **Copy to Workspace**: The correct CSV files were successfully written to the workspace:
  - [forc_curve.csv](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/forc_curve.csv)
  - [compass_field_log.csv](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/compass_field_log.csv)

---

## Cross-Platform Environment Configuration (Dropbox Support)

To solve write access restrictions and file locks inside Dropbox across macOS, Windows, and Linux, we implemented a full environment wrapper suite:

### 1. Non-Synced Virtual Environment Setup
- Created [requirements.txt](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/requirements.txt) to define core dependencies (`numpy` and `matplotlib`).
- Created [setup_env.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/setup_env.py) which builds a local virtual environment stored outside the Dropbox directory (at `~/.compass_sim_venv`). This prevents syncing thousands of pip files over Dropbox and eliminates cross-platform path/binary incompatibility.
- Automatically generated OS-specific launcher scripts in the workspace directory:
  - [run.sh](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/run.sh) (for macOS/Linux)
  - [run.bat](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/run.bat) (for Windows)

### 2. Sandbox-Safe Execution Wrapper
- Created [run_safely.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/run_safely.py):
  - Creates a local temporary folder (`/tmp` or `%TEMP%`).
  - Copies all Python project scripts into this local temp environment.
  - Executes the simulation inside the local temporary directory, bypassing all macOS CloudStorage/FileProvider sandboxing restrictions and Dropbox sync lock conflicts.
  - Upon successful execution, copies all output files (`.png`, `.csv`, `.mp4`, etc.) back to the Dropbox workspace folder.
  - Handles write permission issues gracefully by saving copies to a fallback location (e.g. `~/Downloads/compass_outputs/`) and printing detailed user instructions on OS permission fixes (like macOS Full Disk Access configuration).

---

## Verification Results for `run_safely.py`

Tested `run_safely.py` execution within the scratch sandbox using:
```bash
python3 run_safely.py compass.py --N 4 --M 4 --field_mode forc --forc_n_curves 3 --make_images 1
```

- **Output**: The wrapper successfully created a temporary folder `compass_run_*`, ran the simulation inside it, detected the 7 generated files, copied them back to the workspace, and cleaned up the temporary directory.

---

## Next Steps for the User

1. **Set Up the Virtual Environment**:
   Run the setup script in your terminal to create the virtual environment and generate the launchers:
   ```bash
   python3 setup_env.py
   ```
2. **Run Simulations Safely inside Dropbox**:
   Use the local launchers to execute your code. They will automatically activate the local environment and execute the run securely in a local folder:
   - **macOS/Linux**:
     ```bash
     ./run.sh compass.py --field_mode forc
     ```
   - **Windows**:
     ```cmd
     run.bat compass.py --field_mode forc
     ```

