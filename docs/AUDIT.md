# Codebase Audit — `simula agulhas` (Compass-Needle Lattice Simulator)

**Scope.** This repository is a single-researcher physics-simulation project: 2D lattices of
inertial, dipolar-coupled classical "compass needles," integrated with damped Velocity-Verlet
dynamics, used to study hysteresis, FORC diagrams, avalanches, and demagnetization protocols in
artificial dipolar arrays. It is not a packaged library — there is no `src/` layout, no test
suite, no CI, and the git history is a single "Initial commit" (no real revision history; all
versioning was done by copying whole files).

**Canonical engine.** Of the many simulator files, **`compassV2_2.py`** (1679 lines) is the most
feature-complete and should be treated as canonical going forward. But see §1 and §3.1 — it is
**not** actually the file any of the campaign/sweep scripts invoke, and its own internal metadata
misidentifies it as a different, older file.

**Audit method.** `compassV2_2.py` was read in full directly. A background survey covered the
~150-file `OLD/` archive, the other root-level simulator duplicates, and the 8 secondary
scripts (`avalanche_processor.py`, `calc_defaults_temp.py`, `compass_generate_images.py`,
`plot_default_geometries_clean.py`, `damping_sweep.py`/`V02`/`V03`,
`damping_sweep_analysis.py`, `sweep_damping_hysteresis.py`), cross-checked against
`USER_GUIDE.md`, `USER_GUIDE_compassV02_updated.tex`, `FORC.md`, and `compass_sim_protocol_notes.md`.
No code was changed.

---

## 1. Architecture overview

### 1.1 What actually exists, and which files are live

| File | Role | Status |
|---|---|---|
| `compassV2_2.py` (1679 l) | Most advanced engine: adds "V2.1 hardened" flip/avalanche counters (Schmitt-trigger + dwell) and a numerical-stability guard on top of `compassV02.py`. | **Nominally canonical, but orphaned** — nothing else in the repo imports it (its filename with no `.py`-safe module story aside, nothing references it by name). Only `avalanche_processor.py` depends on the `--event_log` feature it introduced, and does so by reading its **output files**, not by importing it. |
| `compassV2_1 (1).py` (1679 l) | Byte-identical to `compassV2_2.py`. | Dead duplicate-download artifact. |
| `compass.py` / `compassV02_corrected.py` (1458 l, byte-identical) | One generation behind `compassV2_2.py` — has the physics engine, geometry, energy accounting, FORC/demag protocols, but **not** the hardened flip/dwell counters or the `dt_guard` stability monitor. | **This is the file actually imported** (`import compass as cs`) by every sweep/campaign script. |
| `compassV02.py` (1292 l) | Earlier, less-hardened version of the above (179 diff lines from `_corrected`). | Superseded by `compassV02_corrected.py`; kept alongside it for no apparent reason. |
| `compass.py.bkp` / `compass.py.old` (3058 l, byte-identical to each other) | A **different, unrelated, larger lineage** self-identified internally as "V79," with a materially different CLI (`--gpu`, `--video`, `--demag off\|on\|rotational\|anneal`, `--halo_mode`, etc.) and a programmatic API (`make_grid()`, `relax()`, `label_magnetic_domains()`) that was later removed from the `compassV02*`/`compassV2_2` lineage. | Dead, but **this is what `USER_GUIDE.md` documents** — despite the filename implying it's a backup of the current `compass.py`, it is not. |
| `OLD/*.py` (103 files, ~216k lines) | `compass_simT*.py` (Portuguese, LLG-inspired damping — earliest era) → `compass_simV11..V77.py` (English CLI, then GPU/CuPy, then FORC) → `freeze_v58..v77.py` (trivial `shutil.copy2` snapshot scripts). | Fully dead archive; confirmed nothing outside `OLD/` imports it. |

**Net effect:** there are 4 independent lineages of "the simulator" alive in the working tree at
once (V79-orphan, `compassV02*`/`compass.py`, `compassV2_1/2_2`, plus the 103-file `OLD/`
archive), and the one wired up to the automated campaign tooling (`compass.py`) is neither the
most advanced (`compassV2_2.py`) nor the one the user-facing docs describe (`compass.py.bkp`).

### 1.2 Physics model (as implemented in `compassV2_2.py`, and identically in `compass.py`)

Each needle `i` is a point dipole fixed at a 2D lattice site, free to rotate in-plane. Equation of
motion (SI units throughout):

```
I·θ̈ᵢ = τᵢ_dip + τᵢ_ext − b·θ̇ᵢ
```

- **Dipolar torque**: computed from a precomputed pairwise dipole-field tensor (`Axx`, `Axy`,
  `Ayy`, `precompute_dipolar_tensor`, compassV2_2.py:402-456) truncated at a cutoff radius
  (`--cutoff_shells`/`--cutoff_m`), optionally with periodic images (`--pbc`, `--n_images`).
  Each time step is then a matrix-vector product rather than an O(K²) rebuild — the main
  performance-relevant design decision in the codebase.
- **Integrator**: a semi-implicit ("damped") velocity-Verlet scheme
  (`run_simulation`'s per-step block, compassV2_2.py:1237-1244, and the equivalent standalone
  `verlet_step`, compassV2_2.py:1214-1224) that solves the linear damping term implicitly rather
  than explicitly — this is numerically correct and unconditionally stable in the damping
  parameter, not a source of instability. See §3.2 for a **duplication** problem in how it's
  wired up, though.
- **Field protocols**: `static`, `hysteresis`, `forc`, `sine`, `step_up/down`, `pulse`,
  `demag_rot`, `demag_linear` (`FieldProtocol.at`, compassV2_2.py:536-630).
- **Observables**: `Mx`, `My`, `M_proj`, polar order `S1`, nematic/director order `S2`,
  `theta_director`, `q_axis`, dipolar/external/kinetic energies, `omega_rms`/`omega_max`
  (`compute_metrics`, compassV2_2.py:655-708); connected-component "domain" statistics via
  union-find (`domain_statistics`, compassV2_2.py:712-763); two independent avalanche/flip
  counters (drive-axis Schmitt-trigger + rest-angle displacement, both with dwell-time
  debouncing, compassV2_2.py:948-1006 and 1246-1317).
- **Auto-calibration**: if `--moment`/`--inertia` are not given explicitly, they are derived from
  an assumed **rhombus-shaped** steel blade geometry (`compute_moment_from_geometry`,
  `compute_inertia_from_geometry`, compassV2_2.py:169-219 — see §1.1 in the parameter table and
  §4 for a doc mismatch on this point).
- **Outputs**: per-run CSV time series, JSON metadata (physical derived quantities + run config),
  NPZ initial/final states, optional per-event CSV (`--event_log`), and PNG lattice renders.

### 1.3 Secondary tooling built on top

- `avalanche_processor.py` — offline spatio-temporal avalanche clustering (causal union-find
  over `--event_log` output) + power-law/log-normal/exponential MLE fits with a Vuong
  likelihood-ratio test. The only secondary script correctly wired to the canonical engine's
  newest feature (because it only reads output files, not the engine's Python API).
- `damping_sweep.py` → `V02` → `V03` — three generations of a campaign wrapper (hysteresis /
  free-relax / FORC stages across geometry × damping(Q) × seed). **V01 and V02 are
  non-functional** against the current `compass.py` (§3.1, finding B1). Only V03 was migrated
  to the live `run_simulation(args)` API, but it introduces a column-mislabeling bug (§3.1,
  finding B2).
- `damping_sweep_analysis.py` — consumes `damping_sweep*` output: hysteresis-loop metrics
  (area, `Bc`, `Mr`), two avalanche-detection methods, MLE/Vuong statistics, FORC ρ(Br,B)
  diagrams. Written against the V01/V02 CSV convention; silently incompatible with V03 output
  (§3.1, finding B2).
- `sweep_damping_hysteresis.py` — a simpler, single-stage sweep variant. **Non-functional**
  outside `--dry_run`.
- `compass_generate_images.py` / `plot_default_geometries_clean.py` — post-hoc PNG regeneration
  from saved CSV/NPZ state, ~150 lines of drawing code duplicated between the two.
- `calc_defaults_temp.py` — a no-CLI, hardcoded-constants scratch calculator for the
  mass/inertia/moment defaults.

---

## 2. Parameter reference

### 2.1 Physical constants (module-level, `compassV2_2.py:118-131`)

| Constant | Value | Physical meaning | Units |
|---|---|---|---|
| `MU0_OVER_4PI` | `1.0e-7` | μ₀/4π, dipole-field prefactor | T·m/A |
| `STEEL_DENSITY_DEFAULT` | `7850.0` | Needle blade material density | kg/m³ |
| `STEEL_MS_SATURATION_DEFAULT` | `1.59e6` | Saturation magnetization (≈ 2.0 T / μ₀) | A/m |
| `CENTER_DISTANCE_DEFAULT` | `0.013` | Default pivot-to-pivot spacing | m |
| `R_DEFAULT` | `0.0065` (= ½ of the above) | Half-spacing; historical variable kept because `d = 2R` | m |
| `NEEDLE_LEN_DEFAULT` | `0.010` | Default blade length (rhombus long diagonal) | m |
| `NEEDLE_WIDTH_DEFAULT` | `0.003` | Default blade width (rhombus short diagonal) | m |
| `NEEDLE_THICKNESS_DEFAULT` | `0.0004` | Default blade thickness | m |
| `DAMPING_DEFAULT` | `5.0e-8` | Default viscous rotational damping `b` | N·m·s/rad |
| `SOURCE_FILE_TIMESTAMP` | `"2026-07-10T11:50:24-03:00"` | Hardcoded "last updated" string | — |

### 2.2 `compassV2_2.py` CLI parameters (the canonical engine)

**Lattice geometry** (`compassV2_2.py:1581-1589`)

| Param | Meaning | Units | Default | Line |
|---|---|---|---|---|
| `--geometry` | Lattice topology: `square`, `triangular`, `honeycomb` | — | `square` | 1582 |
| `--N` | Rows (or nominal honeycomb height) | count | `16` | 1583 |
| `--M` | Columns (or nominal honeycomb width) | count | `16` | 1584 |
| `--R` | Half the pivot-to-pivot distance (`d = 2R`) | m | `0.0065` | 1585 |
| `--needle_frac` | Legacy blade-length fraction of `2R`; only used with `--use_legacy_size_from_R` | dimensionless | `0.80` | 1586 |
| `--needle_len` | Physical blade length (independent of lattice spacing) | m | `0.010` | 1587 |
| `--needle_width` | Physical blade width | m | `0.003` | 1588 |
| `--use_legacy_size_from_R` | If 1, derive blade size from `needle_frac*2R` instead of explicit `needle_len/width` | 0/1 | `0` | 1589 |

**Needle physical properties** (`1591-1603`)

| Param | Meaning | Units | Default | Line |
|---|---|---|---|---|
| `--moment` | Override magnetic moment `m` (else auto-derived) | A·m² | `None` | 1592 |
| `--inertia` | Override moment of inertia `I` (else auto-derived) | kg·m² | `None` | 1593 |
| `--needle_thickness` | Blade thickness | m | `0.0004` | 1594 |
| `--steel_density` | Blade material density | kg/m³ | `7850.0` | 1595 |
| `--steel_Ms` | Saturation magnetization | A/m | `1.59e6` | 1596 |
| `--steel_Bsat` | Saturation flux density; if set, overrides `steel_Ms` via `Ms = Bsat/μ₀` | T | `None` | 1597 |
| `--pivot_radius` | Pivot hole/cylinder radius | m | `0.0` | 1598 |
| `--pivot_thickness` | Pivot cylinder height | m | `0.0` | 1599 |
| `--pivot_density` | Pivot material density | kg/m³ | `8500.0` | 1600 |
| `--pivot_mass` | Explicit pivot mass override | kg | `None` | 1601 |
| `--damping` | Viscous rotational damping coefficient `b` | N·m·s/rad | `5.0e-8` | 1602 |
| `--damping_noise` | Relative uniform random per-needle damping variation | fraction | `0.0` | 1603 |

**Time integration & avalanche detection** (`1605-1617`)

| Param | Meaning | Units | Default | Line |
|---|---|---|---|---|
| `--t_sim` | Total simulated physical time (ignored for FORC/demag modes, which self-compute duration) | s | `2.0` | 1606 |
| `--dt_factor` | Integration step as a fraction of the natural period `T0` | dimensionless (`Δt/T0`) | `0.04` | 1607 |
| `--noise` | Initial angular noise amplitude (Gaussian) | rad | `1.5` | 1608 |
| `--seed` | RNG seed; time-derived if omitted | int | `None` | 1609 |
| `--log_every` | Write one CSV row every N integration steps | steps | `10` | 1610 |
| `--flip_angle_deg` | Rest-angle displacement threshold for the `flip_angle` channel | deg | `90.0` | 1611 |
| `--flip_band_deg` | Schmitt dead-band half-width around the perpendicular to the drive axis (`flip_field` channel) | deg | `30.0` | 1612 |
| `--flip_dwell_T0` | Dwell time to commit a flip, in units of `T0` | multiples of `T0` | `0.5` | 1613 |
| `--flip_settle_frac` | `|ω|` settling threshold as a fraction of `ω0`, for the `flip_angle` channel | fraction of `ω0` | `0.05` | 1614 |
| `--event_log` | Write a per-committed-event CSV (`step,t,needle_id,channel,theta`) | flag | `False` | 1615 |
| `--dt_guard_alpha` | Stability-monitor threshold on `max_i sqrt(m|Bᵢ|/I)·dt` | dimensionless | `0.35` | 1616 |
| `--dt_guard_substep` | Re-integrate flagged steps with 4 global sub-steps (breaks strict symplecticity) | flag | `False` | 1617 |

**Dipolar cutoff & boundaries** (`1619-1625`)

| Param | Meaning | Units | Default | Line |
|---|---|---|---|---|
| `--cutoff_shells` | Interaction cutoff, in multiples of nearest-neighbor distance `r_nn` | multiples of `r_nn` | `3.5` | 1620 |
| `--cutoff_m` | Absolute interaction cutoff; overrides `cutoff_shells` | m | `None` | 1621 |
| `--pbc` | Enable periodic boundary conditions via finite image sums | flag | `False` | 1622 |
| `--n_images` | Number of periodic images per direction | count | `1` | 1623 |
| `--tensor_mem_limit_gb` | Memory guard on the K×K dipolar tensor | GB | `6.0` | 1624 |
| `--float32` | Use float32 tensor/state (vs. float64) | flag | `False` | 1625 |

**Field protocol** (`1627-1636`)

| Param | Meaning | Units | Default | Line |
|---|---|---|---|---|
| `--field_mode` | `static\|hysteresis\|forc\|sine\|step_up\|step_pos\|step_down\|step_neg\|pulse\|demag_rot\|demag_linear` | — | `static` | 1628 |
| `--B_ext` | Field amplitude; auto-derived from `B_max_factor*B_ref` if omitted | T | `None` | 1629 |
| `--B_max_factor` | `B_ext = factor * B_ref` when `B_ext` omitted | dimensionless | `8.0` | 1630 |
| `--phi_ext_deg` | Field direction relative to +x | deg | `0.0` | 1631 |
| `--field_freq` | Sine-mode frequency | Hz | `1.0` | 1632 |
| `--field_delay` | Delay before step/pulse protocols activate | s | `0.0` | 1633 |
| `--t_pulse` | Pulse duration | s | `None` | 1634 |
| `--hyst_spacing` | `linear` or `log` time-warping of the hysteresis ramp | — | `linear` | 1635 |
| `--hyst_log_k` | Log-spacing concentration parameter | dimensionless | `5.0` | 1636 |

**FORC protocol** (`1638-1644`)

| Param | Meaning | Units | Default | Line |
|---|---|---|---|---|
| `--forc_Br_min` | Minimum reversal field; defaults to `-Bmax` | T | `None` | 1639 |
| `--forc_n_curves` | Number of minor reversal curves | count | `30` | 1640 |
| `--forc_t_sat` | Hold time at +saturation before each curve | s | `0.05` | 1641 |
| `--forc_t_ramp_down` | Ramp-down time to `Br` (used only if `--forc_rate` unset) | s | `0.10` | 1642 |
| `--forc_t_ramp_up` | Ramp-up time back to `+Bmax` (used only if `--forc_rate` unset) | s | `0.20` | 1643 |
| `--forc_rate` | Constant `dB/dt`; overrides fixed ramp times, computing per-curve duration | T/s | `None` | 1644 |

**Demagnetization** (`1646-1649`)

| Param | Meaning | Units | Default | Line |
|---|---|---|---|---|
| `--demag_freq` | Rotation/oscillation frequency of the demag field | Hz | `2.0` | 1647 |
| `--demag_cycles` | Number of decay cycles | count | `20` | 1648 |
| `--t_relax_after` | Extra relaxation time appended after demag completes | s | `2.0` | 1649 |

**Output & performance** (`1651-1662`)

| Param | Meaning | Units | Default | Line |
|---|---|---|---|---|
| `--out_dir` | Output root directory | path | `compassV2_output` | 1652 |
| `--tag` | Run tag (else auto-built from geometry/mode/N/M/seed) | str | `None` | 1653 |
| `--use_gpu` | Use CuPy backend | flag | `False` | 1654 |
| `--progress` | Print progress bar | flag | `False` | 1655 |
| `--verbose` | Print derived-quantity summary | flag | `False` | 1656 |
| `--make_plot` | Generate a quick-look diagnostic PNG | flag | `False` | 1657 |
| `--png_dpi` | Lattice PNG resolution | dpi | `300` | 1658 |
| `--png_transparent` | Transparent PNG background | flag | `False` | 1659 |
| `--png_with_axes` | Include axes/grid/title on lattice PNGs | flag | `False` | 1660 |
| `--png_no_panel_titles` | Suppress panel titles in side-by-side PNG | flag | `False` | 1661 |
| `--domain_tol_deg` | Angular tolerance for domain clustering | deg | `15.0` | 1662 |

**Derived quantities** (computed in `run_simulation`, logged to metadata JSON, not CLI flags)

| Quantity | Formula | Meaning | Units |
|---|---|---|---|
| `B_ref` | `μ₀/4π · 2m / r_nn³` | On-axis dipole field at nearest-neighbor distance — used as the natural field scale | T |
| `B_eff` | `max(|B_ext|, B_ref)` | Field scale used to set the integration timestep | T |
| `omega0` | `sqrt(m·B_eff/I)` | Characteristic (small-oscillation) angular frequency | rad/s |
| `T0` | `2π/omega0` | Characteristic oscillation period | s |
| `dt` | `dt_factor * T0` | Integration timestep | s |
| `Q` | `omega0·I / damping` | Quality factor (damping regime indicator) | dimensionless |
| `u_thresh` | `sin(flip_band_deg)` | Schmitt-trigger threshold on `cos(θ−φ)` | dimensionless |

### 2.3 Secondary-script parameters

Each auxiliary script has its own independent argparse surface (campaign geometry sweeps,
plotting options, analysis thresholds). These are fully enumerated with file:line references in
the companion research notes; the highlights relevant to physical correctness:

| Script | Notable params | Units | Default |
|---|---|---|---|
| `damping_sweep.py`/`V02`/`V03` | `--Q_min`/`--Q_max` (damping quality-factor grid bounds) | dimensionless | `0.05` / `15.0` |
| | `--t_sim_periods` (hysteresis duration) | multiples of `T0` | `40.0` |
| | `--grid_n` (N=M lattice size) | needles/side | `30` |
| `damping_sweep_analysis.py` | `--mad_threshold` (avalanche-detection MAD multiplier) | × MAD | `4.0` |
| `avalanche_processor.py` | `--t_link_T0` (causal clustering time window) | multiples of `T0` | `1,2,4` |
| | `--r_link_rnn` (spatial clustering radius) | multiples of `r_nn` | `1.05` |
| `sweep_damping_hysteresis.py` | `--n_damping` | count | `8` |

---

## 3. Bugs and numerical issues

### 3.1 Correctness bugs (will produce wrong or no results)

**B1 — Two of three campaign scripts call a removed API and cannot run.**
`damping_sweep.py` (lines 165, 171, 235) and `damping_sweepV02.py` (177, 183, 247) call
`cs.make_grid()`, `cs.relax()`, `cs.label_magnetic_domains()` via `import compass as cs`, which
resolves to the root `compass.py`. None of these three functions exist in `compass.py` (or
`compassV02*`/`compassV2_2.py`) — they survive only in the unrelated, dead `compass.py.bkp`
("V79") lineage. **Every stage of both scripts raises `AttributeError` on first use.** Neither
script can currently produce a single data point.

**B2 — The one campaign script that *can* run (`damping_sweepV03.py`) mislabels its own output
column, silently corrupting all downstream avalanche statistics.**
`damping_sweepV03.py`'s `_read_compass_csv_as_field_log()` (~line 207-210) maps the field-log
"S" column to `data["S1"]` — the continuous **polar order parameter** — but
`damping_sweep_analysis.py`'s `detect_avalanches_from_S()` (line ~237-273) assumes "S" is a
sparse, mostly-zero **event count** (`is_active = S > 0`; contiguous nonzero runs = avalanches).
Fed V03 output, this either detects one giant spurious "avalanche" per run or nothing meaningful
at all — it fails silently, producing plausible-looking but physically wrong power-law/Vuong
statistics rather than an error. This is the single highest-priority bug: it can generate
publication-looking numbers that are not measuring what the analysis code claims to measure.

**B3 — `sweep_damping_hysteresis.py` cannot execute at all outside `--dry_run`.**
`run_single()` references module globals `INERTIA`, `B_MAX`, `MOMENT` that are never assigned
anywhere in the file (→ `NameError`); separately, it calls `compute_t_sim(damping, R=R)` against
a signature `def compute_t_sim(R=R_DEFAULT)` — the positional `damping` argument silently binds
to the `R` (needle radius) parameter instead, corrupting the simulated-time calculation even if
the `NameError` were fixed; and it calls `run_single(..., damp_idx=...)` against a signature that
has no `damp_idx` parameter (→ `TypeError`). Any one of these three independently blocks
execution.

**B4 — Metadata self-misidentification: outputs cannot be traced to the code that produced
them.** `compassV2_2.py`'s own docstring header, its JSON metadata (`"program": "compassV02.py"`,
line 1082, `"version": "2.1.0"`, line 1083), and its `SOURCE_FILE_TIMESTAMP` constant (line 131,
hardcoded to `2026-07-10T11:50:24-03:00`, one day *before* the file's actual last edit) all
identify it as `compassV02.py`. Given that four independently-lineaged simulator files coexist in
this repo (§1.1), and the field-log/metadata JSON is the only durable record of which file
produced a given dataset, **every output produced by `compassV2_2.py` records the wrong program
name and a stale timestamp**, making provenance unrecoverable from the data alone.

**B5 — `plot_default_geometries_clean.py:51-52`** (`load_state()`) reads `meta["derived"]`
unconditionally and raises an uncaught `KeyError`/`TypeError` on any `.npz` saved without full
metadata — no fallback, unlike the more defensive loader in `compass_generate_images.py`.

### 3.2 Numerical-design observations (not bugs, but worth flagging)

- **Integrator core is duplicated, not reused.** The per-step update in `run_simulation`'s main
  loop (compassV2_2.py:1237-1244) and the standalone `verlet_step()` helper
  (compassV2_2.py:1214-1224) implement the *exact same* semi-implicit damped-Verlet formula, but
  the main loop never calls `verlet_step()` — it re-derives the same algebra inline, and
  `verlet_step()` is only invoked from the `--dt_guard_substep` re-integration branch. The two
  are currently consistent, but a future correctness fix applied to one (e.g., a damping-noise
  edge case) is not guaranteed to propagate to the other. This is the most safety-relevant
  duplication in the repo because it's the physics kernel itself.
- **The stability guard's correction is off by default.** `--dt_guard_alpha` monitoring is always
  on (free), but sub-stepping correction (`--dt_guard_substep`) is opt-in. In the default
  configuration, steps that violate `omega_local·dt > alpha` are only *counted* and reported as a
  post-hoc warning — the actual trajectory is not corrected unless the user notices the warning
  and reruns. Any run whose warning is not read produces silently under-resolved dynamics with no
  trace in the primary CSV output.
- **Nearest-neighbor "reference field" convention is a physical simplification presented as
  exact.** `B_ref = μ₀/4π · 2m/r_nn³` is the on-axis (aligned) two-dipole field formula, used
  uniformly regardless of lattice geometry or the actual bond-angle distribution (the transverse
  formula would give half this value). It is a reasonable order-of-magnitude normalization, but
  it silently drives the auto-selected `B_ext` (via `--B_max_factor`), the integration timestep
  (`omega0`, `T0`, `dt`), and the "sweep rate in units of `B_ref`/`T0`" reported in metadata —
  all as if `B_ref` were unambiguous, when it is a directional-averaging choice.
- **Dipolar cutoff (default 3.5×`r_nn`) truncates a formally long-range (1/r³) interaction with
  no Ewald-type correction**, and periodic images (`--pbc`, `--n_images`, default 1) are combined
  with the *same* cutoff mask used for the non-periodic case — if `cutoff_m` is smaller than the
  box size, most periodic images are masked out and `--pbc` has near-zero effect; if `cutoff_m`
  exceeds `n_images * L`, real periodic neighbors beyond the first image shell are missed
  entirely. Neither case triggers a warning. For a study whose whole point is long-range dipolar
  correlations and avalanche statistics, an unflagged cutoff/PBC mismatch is a real risk of
  systematic, size-dependent bias.
- **`q_axis` order parameter is computed and logged unconditionally for every geometry**,
  including triangular and honeycomb lattices, even though it measures population imbalance
  between the x/y axes — a concept native to the square lattice's natural easy axes and
  physically ambiguous for a triangular or honeycomb bond geometry. The in-code comment
  acknowledges this ("meaningful for square-type... interpreted cautiously" elsewhere) but the
  CSV column is emitted regardless, inviting misuse by downstream analysis that doesn't know to
  discount it for non-square runs.
- **Fixed global timestep from a two-body neighbor estimate.** `dt` is set once from `omega0`
  (itself from a two-dipole nearest-neighbor field estimate), not from the actual many-body local
  field each needle experiences once the lattice is not in its diluted/aligned reference
  configuration. The `dt_guard` monitor exists precisely to catch this, which is good practice,
  but it means the *default* `dt_factor=0.04` is a heuristic, not a resolved bound — see the
  "guard correction off by default" point above.

---

## 4. Code smells, duplication, and documentation drift

### 4.1 Duplication (see §1.1 for the full inventory)

- **4 parallel simulator lineages** (`compassV2_2`/`compassV2_1(1)`, `compass.py`/`compassV02*`,
  `compass.py.bkp`/`.old` "V79", and the 103-file `OLD/` archive) totaling **~220,000 lines** of
  near-duplicate Python, none reconciled via version control (the repo's single "Initial commit"
  means there is no real diff history — versioning was done entirely by copying whole files,
  including a dedicated `freeze_vNN.py` snapshot-copy convention).
- **~150 lines of PNG-drawing code** (`needle_halves()`, color/style constants, pivot-circle
  drawing) duplicated verbatim between `compass_generate_images.py` and
  `plot_default_geometries_clean.py`.
- **Two near-identical angle-wrapping helpers** in `compassV2_2.py`
  (`wrap_angle`, line 139; `wrap_angle_np`, line 144) plus the same formula inlined a third and
  fourth time in the main integration loop and in `verlet_step()`, instead of all four call sites
  sharing one function.
- **Width/length ratio for the needle blade is defined three different, inconsistent ways** across
  the repo: `0.10` in `USER_GUIDE.md`'s prose description, `0.22` as a fallback guess in
  `compass_generate_images.py:150`, and `0.30` as the actual default/sweep convention
  (`NEEDLE_WIDTH_DEFAULT/NEEDLE_LEN_DEFAULT` and `damping_sweep*.py`'s
  `needle_width = needle_len * 0.30`).
- **`damping_sweepV03.py`'s `_default_compass_args()`** (a hand-maintained namespace shimming the
  engine's ~69-flag schema) does not set any of the six V2.1 flags
  (`flip_band_deg`/`flip_dwell_T0`/`flip_settle_frac`/`dt_guard_alpha`/`dt_guard_substep`/
  `event_log`). If `import compass as cs` were ever repointed at `compassV2_2.py` — the
  nominally-canonical engine — this shim would immediately break with an `AttributeError` on the
  first missing attribute `run_simulation` touches. The campaign tooling and the canonical engine
  are on a collision course by construction.

### 4.2 Documentation drift

- **`USER_GUIDE.md` (the primary onboarding document) documents the dead "V79" `compass.py.bkp`
  lineage, not any live file.** Nearly every CLI flag name differs from `compassV2_2.py`
  (`--gpu` vs `--use_gpu`, `--pbc_images` vs `--n_images`, `--domain_tol` vs `--domain_tol_deg`,
  `--phi_ext` vs `--phi_ext_deg`, entire `--video`/`--frame_every`/`--fps`/`--demag
  off|on|rotational|anneal` subsystems that don't exist in any live engine), and several
  numeric defaults are simply wrong for the current code (`--R` 0.0075 vs actual 0.0065,
  `--dt_factor` 0.05 vs actual 0.04, `--needle_frac` 0.6667 vs actual 0.80, `--pivot_radius`
  0.001 vs actual 0.0). A new contributor following this guide verbatim would construct CLI
  invocations that fail immediately against every live engine file.
- **`USER_GUIDE.md`'s physics description of the needle** ("flat rectangular steel sheet,"
  `I = (1/12)·mass·(L²+W²)`) does not match the actual code, which models a **rhombus**
  (`area = 0.5·L·W`, `I = (1/24)·mass·(L²+W²)`, confirmed by both
  `compute_inertia_from_geometry` and the needle-drawing code in `_needle_halves_for_png`, which
  literally draws two triangles forming a diamond).
- **`USER_GUIDE_compassV02_updated.tex` is well-aligned** with `compassV2_2.py`'s flag set (a
  useful signal that documentation *can* be kept in sync in this project) but predates the V2.1
  hardened-counter era: `--flip_band_deg`, `--flip_dwell_T0`, `--flip_settle_frac`,
  `--dt_guard_alpha`, `--dt_guard_substep`, and `--event_log` (which `avalanche_processor.py`
  depends on) appear nowhere in it.
- **`FORC.md` and `compass_sim_protocol_notes.md`** are, by contrast, accurate design notes — they
  describe the historical "variable sweep rate trap" (V67–V69) and its fix via `--forc_rate`,
  and this narrative matches `compassV2_2.py`'s actual implementation
  (compassV2_2.py:1076-1077). These two documents are not stale; they're just design notes, not
  parameter references, and shouldn't be mistaken for the latter.

### 4.3 Other smells

- Several CLI arguments (`--pivot_radius`, `--pivot_thickness`, `--pivot_density`,
  `--pivot_mass`, `--needle_thickness`, `--steel_density`) have no `help=` text at all in
  `compassV2_2.py`'s `argparse` setup, while adjacent arguments in the same group do — an
  inconsistent documentation bar within one file.
- `calc_defaults_temp.py` has no CLI/argparse surface — every physical input is a hardcoded
  module constant, a workflow outlier relative to every other script in the repo.
- Silent parameter coupling: if geometry (`--needle_len`/`--needle_width`) is changed without
  also setting `--B_ext` explicitly, the auto-derived field scale (`B_ref` → `B_ext`) changes
  nonlinearly and invisibly, which can make parameter sweeps that vary geometry hard to interpret
  without reading the metadata JSON afterward.
- Root directory clutter: dozens of generated artifacts (`.png`, `.mp4`, `.csv`, `.pdf`,
  `.tex`/`.out` LaTeX build products, a 15 MB `results2.zip`) are committed alongside source code
  at the repo root, with no `.gitignore` separating source from generated output (though note:
  the git history shows only one commit, so this may reflect an import rather than an active
  habit — still worth fixing going forward).

---

## 5. Prioritized refactoring plan

### P0 — Fix before trusting any existing sweep-campaign results
1. **Stop using `damping_sweep.py`/`V02`/`sweep_damping_hysteresis.py` entirely** (B1, B3) — they
   cannot run a real simulation today. Do not attempt small patches; they target a removed API.
2. **Quarantine or fix `damping_sweepV03.py` + `damping_sweep_analysis.py` together** (B2) — the
   "S" column semantics must be reconciled (either re-emit a true sparse flip-count column from
   `compass.py`'s CSV, or change `detect_avalanches_from_S` to operate on the real event-log /
   `flip_field` data that `compassV2_2.py` already produces). **Any avalanche statistics already
   produced by this pair of scripts should be treated as unverified until this is resolved.**
3. **Fix the metadata self-misidentification** (B4): update `compassV2_2.py`'s docstring header,
   `"program"`/`"version"` metadata fields, and `SOURCE_FILE_TIMESTAMP` to actually name and date
   itself. Do this for whichever file is chosen as canonical in step 5 below.

### P1 — Consolidate the four simulator lineages into one
4. **Decide, explicitly, which engine is canonical** — today `compassV2_2.py` is the most
   feature-complete but is not what any tooling actually runs; `compass.py`
   (=`compassV02_corrected.py`) is what the campaign scripts import but lacks the hardened
   avalanche counters and stability guard. Recommend promoting `compassV2_2.py` to be the one
   true `compass.py`, then repointing every `import compass as cs` at it and fixing the resulting
   `_default_compass_args()` gaps (the six missing V2.1 flags, §4.1) as part of the same change.
   Do this as a single change, not incrementally, to avoid a period where campaign output is
   silently generated from two different engines.
5. **Archive, don't keep live, the dead lineages**: move `compass.py.bkp`/`.old`,
   `compassV02.py`, `compassV2_1 (1).py`, and the entire `OLD/` directory (521 MB, 103 files) out
   of the active working tree (e.g., into a clearly-labeled `archive/` outside the main project
   folder, or simply delete them given they are fully reconstructable from nothing — there is no
   real git history to lose). This removes ~220k lines of dead code from anyone's search results.
6. **Rewrite `USER_GUIDE.md` against the actual canonical engine's `argparse` output** rather than
   hand-maintained prose — every flag name, default, and physics description in §4.2 needs
   correcting. Consider generating the parameter reference table directly from
   `build_parser()` so it cannot drift again.

### P2 — Reduce duplication that risks future correctness drift
7. **Make the integrator single-sourced**: have `run_simulation`'s main loop call `verlet_step()`
   directly instead of re-deriving the same formula inline, so a future fix only needs to happen
   once.
8. **Collapse `wrap_angle`/`wrap_angle_np`** into one function used at all four call sites.
9. **Extract the needle-drawing code** (`needle_halves`, color constants, pivot circles) shared
   between `compass_generate_images.py` and `plot_default_geometries_clean.py` into one shared
   module.
10. **Resolve the needle width/length ratio** to one number (currently 0.10/0.22/0.30 across
    docs/fallback/actual-default) and state it in exactly one place.

### P3 — Physics/methodology hardening (worth a design discussion, not urgent bugs)
11. Make the `--dt_guard_substep` correction the default whenever the monitor flags any step
    (or, at minimum, fail loudly rather than warn-and-continue when a run is flagged), so
    under-resolved dynamics can't ship unnoticed.
12. Add an explicit warning (or hard validation) when `--pbc` is combined with a
    `cutoff_m`/`n_images`/box-size combination that would make the periodic images a no-op or
    would skip real neighbor shells.
13. Either gate `q_axis` computation/reporting behind `--geometry square`, or add a one-line
    caveat to its CSV column header/metadata note for triangular/honeycomb runs, so downstream
    analysis doesn't average it across geometries where it isn't meaningful.
14. Consider documenting (in one canonical place) the `B_ref` on-axis-vs-transverse convention
    and its downstream effects on auto-selected `B_ext` and reported sweep rates, since it is
    currently presented as an unambiguous physical quantity in metadata JSON.
15. Add a minimal automated smoke test (e.g., `pytest` running each `--field_mode` for a tiny
    4×4 lattice for a handful of steps) — there is currently no test suite at all, so every
    refactor above is unverifiable except by manual inspection.
