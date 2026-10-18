# Implementation Plan — Cross-Platform Dropbox Environment Configuration

This plan details the implementation of a robust, cross-platform environment configuration to run the compass lattice simulation scripts (`compass.py` and sweeps) without any permission or file-locking problems when the project is hosted in a Dropbox (or other cloud-synced) directory on macOS, Linux, or Windows.

## The Problem

1. **macOS CloudStorage/FileProvider Sandboxing**: On macOS, directories under `~/Library/CloudStorage/Dropbox` are managed by FileProvider. Subprocesses (like Python scripts executed inside an IDE terminal) often lack permission to write directly to these folders, raising `PermissionError: [Errno 1] Operation not permitted`.
2. **Cloud Sync File Locking**: On Windows and Linux, Dropbox immediately locks files (plots, CSVs, or frames) the moment they are created or modified to sync them. If the Python script tries to open, write, or rename them at the same time, it crashes with `PermissionError: [Errno 13] Permission denied`.
3. **Dropbox Virtual Environment Syncing**: Storing virtual environments (`.venv`) inside Dropbox is a bad practice. It leads to syncing thousands of tiny dependency files, which degrades system performance, and virtual environments are OS-specific, meaning a synced virtual environment will crash on a different OS.

---

## Proposed Solution

We will create a standardized, automated setup and wrapper script suite:
1. **`setup_env.py`**: A cross-platform setup script that:
   - Detects the user's OS.
   - Installs dependencies (`numpy`, `matplotlib`) into a virtual environment stored **outside** the Dropbox directory (e.g. in `~/.compass_sim_venv/`) so that virtual environment files are never synced.
   - Generates local shell script wrappers (`run.sh` for Mac/Linux and `run.bat` for Windows) to run the simulation using this local venv.
2. **`run_safely.py`**: A robust wrapper that runs Python simulation code in a local temporary directory (e.g., `/tmp` or `%TEMP%`) and copies the output files back. This completely bypasses macOS sandboxing errors and Dropbox file sync locks during simulation execution.
3. **Detailed Setup Guide**: Instructions for system-level settings (such as macOS Full Disk Access).

---

## Proposed Changes

### Configuration & Dependency Setup

#### [NEW] [requirements.txt](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/requirements.txt)
- Define standard packages: `numpy` and `matplotlib`.

#### [NEW] [setup_env.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/setup_env.py)
- Automatically verify Python 3.
- Locate the home directory and create a virtual environment (`.compass_sim_venv`) outside of the Dropbox folder.
- Install packages from `requirements.txt`.
- Generate local scripts in the workspace:
  - `run.sh` (for Mac/Linux)
  - `run.bat` (for Windows)
- Mark `run.sh` as executable (`chmod +x`).

---

### Execution Wrapper

#### [NEW] [run_safely.py](file:///Users/jps/Library/CloudStorage/Dropbox/simula%20agulhas/run_safely.py)
- Accepts the script name and any arguments (e.g., `python3 run_safely.py compass.py --field_mode forc`).
- Creates a temporary directory outside Dropbox.
- Copies the simulation script to the temporary directory.
- Executes the simulation script using `subprocess` in the temporary directory.
- Copies output files (`.png`, `.csv`, `.mp4`) back to the Dropbox folder.
- Handles `PermissionError` during copy operations by offering detailed OS-specific advice (like Full Disk Access instructions) or writing to a fallback location (like the Desktop or Downloads) so that user simulation data is never lost.

---

## Verification Plan

### Manual Verification
1. Run `python3 setup_env.py` to create the environment and ensure the local virtual environment is created outside the Dropbox directory.
2. Verify that `run.sh` (Mac) and `run.bat` (Windows) are generated successfully.
3. Execute the simulation through the wrapper:
   ```bash
   python3 run_safely.py compass.py --N 4 --M 4 --field_mode forc --forc_n_curves 3 --make_images 1
   ```
4. Verify that the outputs are successfully written to the Dropbox workspace directory without any `PermissionError`.
