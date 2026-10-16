# CLAUDE.md — compass-needle lattice simulator

Research code (J. P. Sinnecker, CBPF). The goal is a paper on inertial (Newtonian) dynamics
of 2D arrays of dipolar compass needles. Read `docs/AUDIT_PHYSICS_ROADMAP.md` (physics validity
and roadmap, in Portuguese) and `docs/AUDIT.md` (code organization) before larger changes.

## Layout
- `compass.py` — the only simulation engine. `import compass as cs` everywhere.
- `config.yaml` + `sim_config.py` — all defaults; CLI flags override them.
- `damping_sweepV03.py`, `damping_sweep_analysis.py`, `avalanche_processor.py` — campaign/analysis.
- `tests/` — `test_reference_runs.py` (regression, frozen numbers), `test_hysteresis_protocol.py`,
  `test_physics.py` (physical validity; known defects are `xfail(strict=True)`).
- `archive/` — dead code. Never import from or edit it.
- Work only in this git clone. The Google Drive folder "simula agulhas" is a stale archive.

## Commands
```bash
python3 -m pytest -q tests                      # all tests (~seconds)
python3 compass.py --help                       # CLI
python3 tests/generate_reference.py             # regenerate regression references (deliberately!)
```

## Physics conventions (do not change silently)
- SI units everywhere (T, A·m², kg·m², s). μ0/4π = 1e-7.
- Model: rigid planar rotors, point dipoles at fixed pivots,
  I θ̈ = τ_dip + τ_ext − b θ̇. Integrator: damped velocity Verlet (implicit damping),
  symplectic at b=0. Do NOT replace it with RK4.
- Energy: E_dip = −½ Σ m_i·B_i^dip, E_ext = −Σ m_i·B_ext, E_kin = ½ I Σ ω².
- `B_ref = μ0/4π · 2m/r_nn³` is a normalization, not the local field. The real local field
  in an ordered array is ≈2 B_ref.
- Known defects that Phase 1 fixes (see the roadmap). Keep the matching xfail tests in sync:
  1. Dipolar cutoff 3.5 r_nn biases the field by ~17% and removes edge physics. The tensor is
     dense, so the full sum costs nothing extra. The target default is free boundaries plus the full sum.
  2. dt comes from ω0(B_ref). The local frequency is ~2× higher. Target: dt from max ω_loc.
  3. Q uses B_eff = max(B_ext, B_ref), so it depends on the protocol. Target: an intrinsic Q₀
     from the lattice local field, with the field-dependent Q reported separately.
  4. Sweep rate is not a controlled dimensionless axis. Loop area depends strongly on it.
  5. No work/dissipation ledger. Target: W_field and D = ∫ b ω² dt logged, with the energy balance tested.
- Any change to physics defaults changes the regression references. Regenerate them in the
  same commit and say so in the commit message.

## Scientific framing rules
- This is a macroscopic inertial dipolar-rotor system. It is NOT a model of nanomagnetic artificial spin ice
  (LLG ≠ Newton). ASI is motivation only.
- Never claim avalanche exponents or "universality classes" from these system sizes.
- Any claim about the effect of Q on hysteresis must be made at the quasi-static limit (rate → 0) or
  as an explicit function of the dimensionless sweep rate.

## Workflow
- Branch per phase/fix (`fase1/...`). Small commits. Run `pytest` before committing.
- Campaign outputs go outside the repo (e.g. `~/compass_runs/`). Never commit CSV/MP4/ZIP results.
- Docs for the author are in Portuguese; code, comments and commit messages are in English.
