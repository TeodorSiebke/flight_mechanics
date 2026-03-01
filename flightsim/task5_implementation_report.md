# Task 5: Equations of Motion - Implementation Report

This report documents the implementation of the longitudinal flight simulation model for the Windex aircraft, covering the transfer of aerodynamic data, coordinate transformations, trim analysis, and dynamic exploration.

---

## 1. Aerodynamics Database Transfer
The aerodynamics of the Windex were transferred from the Nonlinear Lifting Line (NLL) model to the flight simulator using a custom generation script.

*   **Database Generation**: `generate_fsim_database.m` was created to sweep Angle of Attack ($\alpha \in [-5^{\circ}, 15^{\circ}]$) and Elevator Deflection ($\delta_e \in [-5^{\circ}, 0, 5^{\circ}]$).
*   **Fixed Reference Point**: For consistency across simulation states, a fixed **Aerodynamic Reference Point (ARP)** was set at **0.25 m** behind the leading edge (LE) of the wing root. This matches the `fsm.pref` datum in the simulator template.
*   **Data Integration**: The results (Tables for $C_L, C_D, C_m$) and stability derivatives ($C_{mq}, C_{l\beta}$, etc.) were hard-coded into `make_fsim.m`.
*   **Files Modified**: `make_fsim.m`, `generate_fsim_database.m`.

## 2. Center of Gravity (CG) Transformation
The simulator was updated to calculate forces and moments about a **Center of Gravity (CG)** that is different from the fixed ARP.

*   **Coordinate mapping**:
    *   **Aerodynamics (Aero)**: $x$-aft positive, $z$-up positive.
    *   **Flight Mechanics (FM)**: $x$-forward positive, $z$-down positive.
*   **Moment Transfer Logic**: Within `fplmod.m`, the pitching moment about the CG ($M_{CG}$) is calculated from the reference moment ($M_{ARP}$) and the force offsets:
    $$\Delta x = x_{CG} - x_{ref} \quad (\text{positive if CG is aft})$$
    $$\Delta z = z_{CG} - z_{ref} \quad (\text{positive if CG is above})$$
    $$M_{CG} = M_{ARP} + \Delta z \cdot X_{force} - \Delta x \cdot Z_{force}$$
*   **Files Modified**: `fplmod.m`, `make_fsim.m`.

## 3. Gliding Trim Implementation
The trim functionality was modified to solve for a steady **unpowered glide** ($\text{Throttle} = 0$).

*   **State Variables**: The solver targets $[u, w, \theta, \delta_e]$ to find equilibrium.
*   **Constraints**:
    *   $\dot{u} = 0, \dot{w} = 0$ (Zero linear acceleration).
    *   $\dot{q} = 0$ (Zero angular acceleration).
    *   $V = V_{set}$ (Specified target airspeed).
*   **Robustness**: A **Damped Newton-Raphson** loop was implemented in `main.m` to ensure convergence. A damping factor (0.5) and step limits were added to prevent divergence in high-alpha regimes.
*   **Files Modified**: `main.m`, `main_tutorial.m`.

## 4. Glide Polar Analysis
The updated trim function was used to automate a sweep of airspeeds from 23 m/s to 45 m/s.

*   **Behavior**: The simulation captures the transition from high Alpha/high sink (near stall) to minimum sink (approx 27 m/s) to high-speed glide.
*   **Comparison**: The resulting black polar curve in the simulator matches the physical trends of the previous NLL work but accounts for 6-DOF gravitational and mass interactions.
*   **Stall Limit**: The simulation correctly identifies the inability to trim below ~21 m/s (approx 12° Alpha).

## 5. Time Domain Exploration (ODE45)
A dynamic exploration script, `explore_controls.m`, was created to analyze the aircraft's stability using `ode45`.

*   **Control Implementation**: `fsm.de_set` defines the timeline of elevator inputs.
*   **Observed Modes**:
    *   **Short Period**: Rapid damping of pitch oscillations following a "doublet" input, confirming high tail efficiency ($C_{mq}$).
    *   **Phugoid**: Long-period oscillations (energy exchange between speed and altitude) following a pulse input. Due to high $L/D$, these oscillations take >100s to dampen, which is characteristic of sailplanes.
*   **Files Created**: `explore_controls.m`, `main_tutorial.m`.

---

### **Summary of Script Functionality**
-   **`make_fsim.m`**: The "Blueprint" (Data storage).
-   **`fplmod.m`**: The "Physics Engine" (Force/Moment calculations).
-   **`simsub.m`**: The "Pilot" (Control timeline during ODE integration).
-   **`main.m` / `explore_controls.m`**: The "Simulation Control" (Running the solver).
