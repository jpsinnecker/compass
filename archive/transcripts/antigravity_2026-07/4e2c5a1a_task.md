# FORC Integration Tasks

## damping_sweepV02.py
- [x] Add `run_one_forc()` function
- [x] Add `save_forc_csv()` function
- [x] Extend `run_job()` with stage='forc' branch
- [x] Add `forc_dir` output directory in main()
- [x] Add FORC CLI flags
- [x] Add forc jobs to job list
- [x] Update quick_test override (forc_n_curves=3)
- [x] Update MANIFEST_FIELDS

## damping_sweep.py (serial V01)
- [x] Mirror all V02 FORC additions in serial form

## damping_sweep_analysis.py
- [x] Add `load_forc_data()`
- [x] Add `compute_forc_distribution()`
- [x] Add `plot_forc_diagram()`
- [x] Add `plot_forc_vs_Q()`
- [x] Extend main() FORC stage loop + progress bar
- [x] Add forc_diagrams.png to output list
- [x] Remove --forc_dir stub

## Verification
- [ ] Quick-test pipeline end-to-end
- [ ] Confirm forc_diagrams.png produced
