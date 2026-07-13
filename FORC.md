# FORC Field Sweep Rate Protocol for `compass_sim`

This document outlines the theoretical reasoning and practical implementation strategy for defining a physically rigorous magnetic field sweep rate ($dH/dt$) within the `compass_sim` framework. 

---

## 1. The Physics Context & Challenge

Unlike standard quasi-static or Monte Carlo simulations where time is an implicit or arbitrary iteration counter, `compass_sim` implements **real Newtonian inertial dynamics** in SI units:

$$I\ddot{\theta}_{i} = \tau_{i}^{\text{dip}} + \tau_{i}^{\text{ext}} - b\dot{\theta}_{i}$$

Because physical time ($t_{\text{sim}}$) is solved explicitly using a Velocity-Verlet integration scheme, the rate at which the external field changes ($dB/dt$) directly governs the non-equilibrium relaxation pathways of the system. 

If the field is swept too quickly:
- Magnetic moments (needles) cannot overcome local energy barriers in time, causing an artificial inflation of the coercive field ($B_c$).
- In high-$Q$ (underdamped) regimes ($Q \gg 1$), rapid sweeps will induce heavy mechanical oscillations that mask the true avalanche scaling behavior.

---

## 2. Analysis of the Current Framework

As of version V67–V69, the First-Order Reversal Curve (FORC) implementation (`--field_mode forc`) defines the timeline of each minor loop using explicit time intervals rather than a direct rate:
- `B_ext` ($B_{\text{max}}$): The absolute saturation field.
- `--forc_Br_min` ($B_{r,\text{min}}$): The lower bound for the reversal fields.
- `--forc_t_ramp_down`: The time allocated to ramp down from $+B_{\text{max}}$ to a given reversal field $B_{r,k}$.
- `--forc_t_ramp_up`: The time allocated to sweep from $B_{r,k}$ back up to $+B_{\text{max}}$.

### The "Variable Sweep Rate" Trap
If `--forc_t_ramp_up` ($t_{\text{up}}$) is held constant for every reversal curve across the suite, the average sweep rate ($R_k$) for the $k$-th curve becomes a function of its depth:

$$R_k = \frac{dB}{dt} = \frac{B_{\text{max}} - B_{r,k}}{t_{\text{up}}}$$

Consequently, curves with deep reversal fields ($B_{r,k} \approx B_{r,\text{min}}$) will be swept **drastically faster** than minor loops resting near positive saturation. For dynamic systems with inertia, this introduces an uncontrolled variable into your avalanche exponent ($\tau$) calculations.

---

## 3. Recommended Protocol for a Constant Sweep Rate

To ensure rigorous statistical mechanics comparisons (especially when evaluating how the quality factor $Q$ alters the avalanche universality class), the sweep rate must remain perfectly constant across all curves. Two methods are available to enforce this:

### Method A: Orchestration Script Wrapper (No Core Code Modification)
If you manage the simulation programmatically via a script like `damping_sweepV03.py` (the campaign wrapper actually compatible with the current engine — `damping_sweep.py`/`V02` call a removed API and cannot run at all, see `docs/AUDIT.md` bug B1), you can dynamically scale the parameter `--forc_t_ramp_up` for each curve cycle $k$. In practice this is exactly what `damping_sweepV03.py`'s own `--forc_rate` flag already does.

For a target constant sweep rate $R$ (in Tesla/second), calculate the required time parameter for each curve instance before calling the simulator:

$$t_{\text{ramp\_up}, k} = \frac{B_{\text{max}} - B_{r,k}}{R}$$

---

### Method B: True Constant $dB/dt$ Core Modification
Modify the internal field update loop inside `compass.py` for `field_mode == "forc"`. Instead of computing time increments from a fixed total duration, specify a target sweep rate $R$ and update the field linearly step-by-step:

```python
# Conceptual implementation inside the Verlet loop for the up-branch:
B_ext_current += target_sweep_rate * delta_t
```

This decouples the simulation path from predefined segment times, keeping $dB/dt \equiv R$ perfectly invariant.

---

## 4. Physical Calibration & Constraints

When configuring your target constant sweep rate $R$, ensure it satisfies the following dimensionless boundaries relative to your system geometry:

1. **Adiabaticity relative to the Natural Period ($T_0$):**
   The field should change negligibly over a single natural oscillation period of the needle.

   $$R \cdot T_0 \ll B_{\text{ref}}$$

   Where $T_0 = 2\pi / \omega_0$ and $B_{\text{ref}}$ is the reference nearest-neighbor dipolar field.

2. **Damping ($Q$) Coupling:**
   In underdamped regimes ($Q = 15$), a continuous field sweep acts as a periodic driving force that can resonate with the lattice. Ensure the sweep rate is slow enough to allow local torques to settle, avoiding purely inertial cascading artifacts that do not represent true magnetic field-driven avalanches.

---

## 5. The Same Principle, Applied to Hysteresis Mode

The variable-sweep-rate trap from Section 2 is not FORC-specific: any protocol that derives $dB/dt$ implicitly from a fixed total duration over a variable-length path is at risk of it.

`--field_mode hysteresis` has since gained its own rate-varying feature, `--hyst_slow_window B_lo,B_hi` + `--hyst_slow_factor f`, which divides the sweep rate by $f$ while $|B|$ is inside $[B_{\text{lo}}, B_{\text{hi}}]$ (both branches). This is implemented following Method B above, not Method A: `build_hysteresis_schedule()` constructs an explicit piecewise-linear `HystSchedule` of `HystLeg`s (constant rate per leg, split at the exact window boundaries), and the branch/global durations are recomputed from those legs rather than held fixed. Concretely:

- Each leg's rate is a first-class, explicit quantity (`HystLeg.rate_T_per_s`), never inferred after the fact from a fixed total time divided across unequal segments.
- The total simulated time (`t_sim` used for `n_steps`) is derived from the sum of actual leg durations, so it grows to accommodate a slow window rather than silently compressing the rate outside it.
- Every run's metadata JSON records the full segment list (`t_start`, `t_end`, `B_start`, `B_end`, `rate_T_per_s`, `branch` per leg) plus the base and slowed rates, so — exactly as Section 2 argues for FORC curves — the actual rate a dataset was produced under is always explicit and reported, never left to be reverse-engineered from `t_sim / n_segments`.

With no window requested (the default), the hysteresis schedule collapses to the original five-equal-duration branches at one constant rate, so this generalization changes nothing for existing datasets.