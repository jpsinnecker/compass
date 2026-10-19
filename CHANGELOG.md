# Changelog

All notable changes to the simulator (`compass.py`) and to the scientific report.
Format inspired by [Keep a Changelog](https://keepachangelog.com/). Full history,
including the pre-git period (15/06–11/07/2026), in `docs/CRONOLOGIA.md` and on
the branch `history/pre-git`.

## Versioning policy

* **Engine:** semantic versioning in `compass.py` (`__version__`), recorded in every run's metadata (`version`).
  * **MAJOR:** incompatible change of the model or of the output format.
  * **MINOR:** new protocol or observable, or a change of physics defaults (results change).
  * **PATCH:** fixes that do not change the numbers of a correct run.

  Bump the version in the same commit as the change, and tag the commit `vX.Y.Z` once it is on `main`.
* **Scientific report:** `Report Vnn` (`docs/project_report.tex`), numbered continuously since July (V01–V53, then V60, V61). The next revision is V62. Tag `report-Vnn`.
* **Pre-git versions** keep their historical names (T1–T10, V11–V79, V2.0, V2.0.2, V2.1), with tags `code-T1` … `code-V2.1` on `history/pre-git`.
* **Do not confuse** "code V60" (05/07/2026, CLI translated to English) with "Report V60" (17/10/2026).

## [Unreleased] (branch `fase5/bench-protocol`)
### Added
- Bench protocol for a single real needle (`docs/BANCADA.md`, in Portuguese) and its analysis script `bench_fit.py`: I from geometry or added inertia, m/I from frequency versus field, m from the dipole field, viscous damping and Coulomb friction from a fit of the full equation of motion to a free decay, static friction from the dead band. Tests with synthetic data in `tests/test_bench_fit.py`. The engine is unchanged.

## [2.4.0] – 2026-10-18 (branch `fase2/dynamic-hysteresis`, PR #7)
### Added
- Event-driven quasi-static loop `--field_mode hysteresis_qs`: field steps of `--qs_dB_factor`·B_ref, relaxation at fixed field until rest, one CSV row per step (per-step avalanche size and dissipation). The field jump is booked exactly in the energy ledger.
- `--qs_dynamics overdamped`: gradient flow b θ̇ = τ (Q₀ → 0 limit).
- **Saddle check** `--qs_kick` (default 1e-4 rad): after the first rest, kick and relax again before logging; kick energy booked as `E_kick`.
- `--qs_state_log` and `--qs_stop_after`, for targeted reruns.
- Campaign and analysis tools: `fase2_campaign.py` (sweep/QS grids, sharding across machines, GPU mode), `fase2_analysis.py`, `fase2_qs_analysis.py`, `qs_path_selection.py`.
### Fixed
- The QS rest criterion accepted zero-torque saddles, so the gradient flow stayed on unstable symmetric branches. The first QS grid (made without the kick) is invalid for comparing dynamics (`docs/QS_SELECAO_DE_CAMINHO.md`).
### Docs
- User manual with a single Markdown source and a pandoc-built PDF; Report V60 restructured, Report V61 with Results III from the corrected quasi-static grid; `docs/CRONOLOGIA.md`.

## [2.3.0] – 2026-10-16 (Phase 1; PRs #4–#6, on `main`; tag `phase1-done`)
### Changed (results change)
- Full dipolar pair sum by default (`--cutoff_shells 0`). The old 3.5 r_nn cutoff biased the interior field by 15–17 %.
- Time step from a worst-case local-field bound: `dt_factor` is now relative to T_loc, and ω_loc·dt ≤ 0.157.
- Metadata `Q` is the intrinsic Q₀ = √(m B_lat I)/b; `Q_ref` and `Q_field` are kept for comparison.
### Added
- Energy ledger: CSV columns `W_field_J` and `E_diss_J`, plus the metadata block `energy_ledger`.
- Dimensionless sweep rate `--sweep_rate_dimless` (r = (dB/dt)·τ_relax/B_ref), recorded for every hysteresis/FORC run.
- Physics-validity tests (`tests/test_physics.py`), README, `requirements.txt`, MIT license.

## [2.2] – 2026-07-12 (tag `v2.2`, set in July)
- Hardened counters and dt guard as the working engine (`compassV2_2.py`, byte-identical to V2.1); repository consolidation.

## [2.1.0] – 2026-07-12 (tag `code-v2.1.0-canonical`)
- V2.1 promoted to the canonical `compass.py`.
- `config.yaml` with a typed loader, regression tests.
- Hysteresis slow window and adaptive logging (`--hyst_slow_window`, `--hyst_slow_factor`, `--log_adaptive`).

## Pre-git (15/06–11/07/2026)
Summarized in `docs/CRONOLOGIA.md`. Every archived version is a commit on `history/pre-git`:
- **V2.1** (11/07): ringing-proof flip counters, stability monitor, `--event_log`.
- **V2.0.2 / V2.0** (10/07): research branch rewritten from scratch.
- **V79** (09/07): physical needle geometry (0.26 mm, pivot).
- **V78** (08/07): M HUD box, time-based seed, 10 s countdown (recovered from the editor's local history).
- **V67–V77** (06–07/07): FORC and demagnetization.
- **V55–V58** (05/07): V57 reconstructed (never saved; the freeze failed) and the original V55/V56 restored (the archived files were edited in place), in `archive/OLD/reconstructed/`. Snapshots from V58 on were made by `freeze_vNN.py` helpers.
- **V40** (04/07): precomputed dipolar tensor.
- **Campaign scripts** (≤ 05/07): earliest surviving `damping_sweep*.py`, `damping_sweep_analysis.py`, `sweep_damping_hysteresis.py` in `archive/OLD/campaign_early/`.
- **V20** (17/06): formal version record begins.
- **T1** (15/06): first script.

## Scientific report
| Version | Date | Note |
|---|---|---|
| V01–V32 | 04–05/07 | `relatorio_projeto_V01…V22` files (Portuguese, then English from V12); V01–V09 recovered from `OLD/ALLfiles` (tags `report-V01-pt` … `report-V09-pt`) |
| V33–V50 | 06–09/07 | `project_report.tex` in the Antigravity sessions (documents code V66–V79); only V47 survived (tag `report-V47`) |
| V51, V52 | 10/07 | rewritten around the V2 branch (tags `report-V51`, `report-V52` on `history/pre-git`) |
| V53 | 12/07 | git, canonical engine (tag `report-V53`) |
| V60 | 17/10 | restructured as a methods-and-results report (tag `report-V60`); its Results III contained saddle artefacts |
| V61 | 18/10 | Results III, abstract and discussion rewritten with the corrected quasi-static grid (tag `report-V61`) |
