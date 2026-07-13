# FORC and Demagnetization Protocols for Magnetic Compass Arrays

This document aggregates the physical reasoning, mathematical framework, and programmatic implementation strategies discussed regarding First-Order Reversal Curve (FORC) and demagnetization protocols for the `compass_sim` physical engine.

---

## Part 1: Defining Field Sweep Rates in FORC Protocols

In First-Order Reversal Curve (FORC) measurements, defining the magnetic field sweep rate ($dH/dt$) is a critical balancing act between technical limitations (instrumental lag, eddy currents) and physical relaxation mechanisms (magnetic viscosity).

### 1. Mathematical Framework
Depending on your framework, the sweep rate can be implemented in one of two modes:
* **Continuous Sweep Mode:** The applied field is ramped continuously at a fixed rate:
$$\frac{dH}{dt} = R$$
* **Step-by-Step Mode (Settle Mode):** The field is modified in discrete steps ($\Delta H$) and held for a stabilization pause time ($\tau$), giving an effective sweep rate of:
$$R_{\text{eff}} = \frac{\Delta H}{\tau}$$

### 2. Python Implementations
* **Dynamic/Kinetic Approach (e.g., LLG or Kinetic Monte Carlo):** When physical time ($t$) exists explicitly, the field increment per time-step ($\Delta t$) must match the target rate:
$$\Delta H = R \cdot \Delta t$$
* **Quasi-Static Approach (e.g., Stoner-Wohlfarth):** If the system relaxes instantly to a local energy minimum without explicit time dependence, the "sweep rate" is strictly governed by the density of the field grid arrays ($\Delta H$).

---

## Part 2: Implementation & Calibration for `compass_sim`

### 1. The Physics Context
The `compass_sim` simulator solves real Newtonian inertial rotational dynamics in SI units [cite: 6]:

$$I\ddot{\theta}_{i} = \tau_{i}^{\text{dip}} + \tau_{i}^{\text{ext}} - b\dot{\theta}_{i}$$

Because physical time ($t_{\text{sim}}$) is resolved step-by-step using a Velocity-Verlet integrator [cite: 6, 34], the field sweep rate ($dB/dt$) acts as a physical control parameter that can shift the system's dynamic relaxation pathways and change avalanche statistics [cite: 9, 20].

### 2. Resolving the "Variable Sweep Rate Trap"
In older configurations (V67–V69), minor reversal loop timing was defined via fixed segment durations (`--forc_t_ramp_up`) [cite: 160]. This meant the average sweep rate for the $k$-th minor loop depended on its reversal field depth ($B_{r,k}$) [cite: 447, 449]:

$$R_k = \frac{B_{\text{max}} - B_{r,k}}{t_{\text{ramp\_up}}}$$

This approach sweeps deep curves drastically faster than shallow ones. To solve this, version V71 introduced the `--forc_rate` ($R$) parameter [cite: 163, 349], which pre-computes a variable time duration for each unique branch to enforce a perfectly flat field sweep rate across the entire suite [cite: 163, 449]:

$$t_{\text{ramp\_up}, k} = \frac{B_{\text{max}} - B_{r,k}}{R}$$

### 3. Calibration Constraints
To ensure a clean scientific evaluation, your selected sweep rate $R$ must respect two physical boundaries:
1. **Adiabaticity Limit:** The field must change negligibly during one natural oscillation period ($T_0$) of a needle within the reference neighbor dipolar field ($B_{\text{ref}}$) [cite: 35, 37]:
$$R \cdot T_0 \ll B_{\text{ref}}$$
2. **Damping ($Q$) Coupling:** In underdamped regimes ($Q \gg 1$), sweeping the external field too quickly acts as an unwanted mechanical driving force, introducing resonant oscillations that mask the true field-driven avalanche exponents ($\tau$) [cite: 41, 43, 203].

---

## Part 3: Evaluating & Improving the Demagnetization Protocol

> **Note (added after the fact):** the analysis below describes `compass.py` "V70"'s
> `--demag on` alternating-field option, which no longer exists in the codebase — it
> belonged to an older lineage, since archived (see `docs/AUDIT.md`, `USER_GUIDE.md`).
> Of the two methods proposed below, **only Method A was adopted**: its formula is
> exactly what `--field_mode demag_rot` implements today. **Method B was not adopted** —
> the current `--field_mode demag_linear` is a simpler linear-decay variant along a fixed
> direction, with no thermal-noise/annealing mechanism; `compass.py` has no stochastic
> torque term at all.

### 1. Analysis of the Current Alternating Field Demagnetization (AFD)
The framework in `compass.py` (V70) uses standard linear Alternating Field Demagnetization via the `--demag on` option [cite: 162, 453]. The field oscillates along a fixed, single direction ($\phi_{\text{ext}}$) while its peak envelope drops aggressively by 30% per cycle [cite: 162, 454]. This method faces two major physical drawbacks in a 2D interacting dipole grid:
* **Directional Anisotropy:** Because the shaking force is strictly uniaxial, components of magnetic needles perpendicular to $\phi_{	ext{ext}}$ do not experience full resetting torques, leaving a directional texture behind.
* **Coarse Geometric Decay:** A 30% reduction per cycle drops too fast. Highly correlated long-range dipolar lattices easily get trapped in unwanted metastable energy valleys, leaving behind large residual domains instead of randomizing completely.

### 2. Suggested Methods for a Truer Isotropic Ground State

#### Method A: Biaxial Rotational Demagnetization (Recommended)
Instead of shaking along a single line, the external field vector rotates continuously in the $xy$-plane while its magnitude shrinks continuously to zero.

* **Mathematical Form:**
$$B_x(t) = B_{\text{max}} \cdot \left(1.0 - \frac{t}{t_{\text{demag}}}\right) \cdot \cos(2\pi \cdot f_{\text{demag}} \cdot t)$$
$$B_y(t) = B_{\text{max}} \cdot \left(1.0 - \frac{t}{t_{\text{demag}}}\right) \cdot \sin(2\pi \cdot f_{\text{demag}} \cdot t)$$
* **Advantage:** Dragging the field vector spherically or circularly through the plane forces all needles through a complete angular sweep, erasing any directional bias and providing a truer randomized state ($S \rightarrow 0$).

#### Method B: Linear Sweep with Simulated Thermal Annealing
Leverage the mechanical engine's inertial nature by introducing a temporary random thermal noise torque (white noise) during a slow linear field ramp down, simulating a high material temperature ($T$). Ramping down this artificial thermal noise alongside the field allows configurations to bounce out of shallow local minima and reach an isotropic global ground state.
