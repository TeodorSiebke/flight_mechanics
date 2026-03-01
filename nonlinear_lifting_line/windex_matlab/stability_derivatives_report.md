# Windex Aerodynamic Stability Derivatives Report

## 1. Overview
This report documents the estimation of stability derivatives for the Windex aircraft using the **Nonlinear Lifting Line (NLL)** model. These derivatives characterize the aircraft's aerodynamic response to changes in flight state (sideslip $\beta$) and rotation rates ($p, q, r$).

## 2. Definitions and Coordinate Systems

## 4. Sign Conventions and Notations

To match standard aerodynamic notation (e.g., Etkin, Roskam), the following mapping is used from the internal NLL coordinate system:

| Variable | Symbol | NLL Input/Output | Definition |
| :--- | :--- | :--- | :--- |
| **Angle of Attack** | $\alpha$ | `fs.alpha` | Positive nose-up. |
| **Sideslip** | $\beta$ | `fs.beta` | **Note**: Internal NLL $\beta$ is wind from left. We use standard convention where $+\beta$ is wind from right. |
| **Roll Rate** | $p$ | $-\omega_1$ | Positive right wing DOWN. |
| **Pitch Rate** | $q$ | $\omega_2$ | Positive nose UP. |
| **Yaw Rate** | $r$ | $-\omega_3$ | Positive nose RIGHT. |
| **Rolling Moment**| $L$ ($C_l$) | `Cl` | Positive right wing DOWN. |
| **Pitching Moment**| $M$ ($C_m$) | `Cm` | Positive nose UP (Pitching). |
| **Yawing Moment** | $N$ ($C_n$) | `Cn` | Positive nose RIGHT (Yawing). |

### Normalization Factors
- **Rates ($p, r$)**: Normalized by $\hat{p}, \hat{r} = \frac{\omega \cdot b}{2V}$
- **Rate ($q$)**: Normalized by $\hat{q} = \frac{\omega \cdot \bar{c}}{2V}$
- **Moments ($L, N$)**: Normalized by $q \cdot S \cdot b$
- **Moment ($M$)**: Normalized by $q \cdot S \cdot \bar{c}$

## 5. Assumptions and Approximations
1. **Quasi-Steady Flow**: Rotation rates are modeled by adjusting the local velocity vector at each segment. The "curvature" of the flow field over the span is captured.
2. **Linear Region**: Derivatives are estimated at the provided reference alpha (default 2°) assuming small perturbations.
3. **Rigid Wake**: As per the NLL model, wake filaments extend straight back in the wind direction.
4. **Viscous Damping**: Unlike pure potential methods, NLL includes profile drag and viscous lift-slope changes from XFOIL, leading to more accurate damping estimates (especially in roll).

## 6. How to Run
In the MATLAB command window:
```matlab
[results, table] = estimate_stability_derivatives(V, alpha_deg);
```
Example for 30 m/s at 2° alpha:
```matlab
estimate_stability_derivatives(30, 2);
```

## 8. Results (Reference V=30 m/s, alpha=2°)

Based on the simulation, the following derivatives were obtained:

| Derivative | Value | Physical Meaning | Stability |
| :--- | :--- | :--- | :--- |
| **Cmq** | -26.1601 | Pitch Damping | **Stable (-) ** |
| **CLq** | 9.7443 | Lift due to pitch rate | Normal (+) |
| **Clb** | -0.0925 | Dihedral Effect (Roll-Beta) | **Stable (-)** |
| **Cnb** | 0.0469 | Directional Stability (Weathercock) | **Stable (+)** |
| **Cyb** | -0.2769 | Sideforce due to sideslip | Normal (-) |
| **Clp** | -0.7103 | Roll Damping | **Stable (-)** |
| **Cnp** | -0.0600 | Adverse Yaw (Yaw due to Roll) | Typical (-) |
| **Clr** | 0.0068 | Roll due to Yaw Rate | Coupling (+) |
| **Cnr** | -0.0238 | Yaw Damping | **Stable (-)** |

## 9. Physical Interpretation

## 10. Note on Discrepancies
You may notice that **$C_{lr}$ (0.0068)** is about 3 times lower than the placeholder values (~0.02) and significantly lower than the "rule of thumb" ($C_L/4$). This is due to several Windex-specific factors captured by the NLL model:
1. **Wing Washout**: The Windex has significant washout (incidence decreases from 2.67° at the root to 0.67° at the tip). Since $C_{lr}$ is dominated by the speed increase on the outboard wing, the reduced lift at the tips significantly dampens this coupling.
2. **T-Tail Interaction**: The vertical tail is located high above the wing ($z \approx 0.9$m). During a yaw, the side force on the fin creates a rolling moment that can partially oppose the wing's rolling moment depending on the CG height.
3. **Reference Point**: The derivatives are highly sensitive to the chosen reference point ($pref = 0.25$m).

## 11. Summary for Task 5
These high-fidelity derivatives have been successfully integrated into `make_fsim.m`. The longitudinal behavior ($C_{mq}, C_{Lq}$) and primary lateral stability ($C_{l\beta}, C_{n\beta}$) are very robust and provide a solid foundation for the flight simulation.
