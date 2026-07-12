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

As of version V67–V69, the First-Order Reversal Curve (FORC) implementation (`-field_mode forc`) defines the timeline of each minor loop using explicit time intervals rather than a direct rate:
- `B_ext` ($B_{\text{max}}$): The absolute saturation field.
- `-forc_Br_min` ($B_{r,\text{min}}$): The lower bound for the reversal fields.
- `-forc_t_ramp_down`: The time allocated to ramp down from $+B_{\text{max}}$ to a given reversal field $B_{r,k}$.
- `-forc_t_ramp_up`: The time allocated to sweep from $B_{r,k}$ back up to $+B_{\text{max}}$.

### The "Variable Sweep Rate" Trap
If `-forc_t_ramp_up` ($t_{\text{up}}$) is held constant for every reversal curve across the suite, the average sweep rate ($R_k$) for the $k$-th curve becomes a function of its depth:

$$R_k = \frac{dB}{dt} = \frac{B_{\text{max}} - B_{r,k}}{t_{\text{up}}}$$

Consequently, curves with deep reversal fields ($B_{r,k} \approx B_{r,\text{min}}$) will be swept **drastically faster** than minor loops resting near positive saturation. For dynamic systems with inertia, this introduces an uncontrolled variable into your avalanche exponent ($\tau$) calculations.

---

## 3. Recommended Protocol for a Constant Sweep Rate

To ensure rigorous statistical mechanics comparisons (especially when evaluating how the quality factor $Q$ alters the avalanche universality class), the sweep rate must remain perfectly constant across all curves. Two methods are available to enforce this:

### Method A: Orchestration Script Wrapper (No Core Code Modification)
If you manage the simulation programmatically via a script like `damping_sweep.py`, you can dynamically scale the parameter `-forc_t_ramp_up` for each curve cycle $k$.

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