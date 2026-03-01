# Windex Surface Analysis Documentation

This document summarizes the methodology and results for the individual aerodynamic analysis of the Windex lifting surfaces (Wing and Horizontal Tail) using the Nonlinear Lifting Line Theory (NLLLT) solver.

## Analysis Methodology

To find the precise Aerodynamic Center (AC) of each surface, we implemented two dedicated scripts: `analyze_wing.m` and `analyze_tail.m`. These scripts follow a rigorous numerical approach to ensure compatibility with standard aerodynamic theory.

### 1. Geometric MAC Integration
Rather than using a simplified trapezoidal approximation, we numerically integrate the planform properties directly from the NLL geometry definition.
- **Area ($S$)**: $\sum (c_i \cdot |dy_i|)$
- **MAC ($c$)**: $\sum (c_i^2 \cdot |dy_i|) / S$
- **X LE MAC**: The longitudinal position of the MAC Leading Edge is found by weighted average of the segment leading edges.

### 2. Vertical Reference Alignment (Z-Correction)
A key discovery was the **Z-couple effect**. If the moment reference point is at $z=0$ but the lift force acts at a vertical offset (e.g., $z = 0.98$ for the tail), the slight backward tilt of the lift vector at higher $\alpha$ creates a significant pitching moment ($M_{z-couple} \approx -z \cdot L \sin\alpha$). 

To isolate the purely geometric AC, we aligned the reference point $z_{ref}$ with the **Mean Z-height** of each surface.

### 3. Linear Slope Analysis
We Perform an Alpha sweep and extract $C_L$ and $C_m$ (referenced to the MAC LE). We fit a linear regression in the **-2° to 5°** range to find the slopes:
- $C_{L\alpha} = \frac{dC_L}{d\alpha}$
- $C_{m\alpha} = \frac{dC_m}{d\alpha}$

The Aerodynamic Center (Neutral Point) is then calculated as:
$$h_{ac} = - \frac{C_{m\alpha}}{C_{L\alpha}}$$

## Results Summary

Table 1: Geometric Properties
| Surface | Area ($S$) | MAC ($c$) | $X_{LE,MAC}$ | $Z_{Ref}$ (Mean) |
| :--- | :--- | :--- | :--- | :--- |
| **Main Wing** | 7.2537 m² | 0.6274 m | 0.0468 m | 0.1411 m |
| **H-Tail** | 0.8600 m² | 0.4149 m | 2.8166 m | 0.9800 m |

Table 2: Aerodynamic Centers (Ref @ MAC LE)
| Surface | $C_{L\alpha}$ (/rad) | $C_{m\alpha}$ (/rad) | $dC_m/dC_L$ | **AC (% MAC)** |
| :--- | :--- | :--- | :--- | :--- |
| **Main Wing** | 6.0411 | -1.6532 | -0.2737 | **27.37%** |
| **H-Tail** | 4.9535 | -1.3195 | -0.2664 | **26.64%** |

> [!NOTE]
> The deviation from the theoretical 25% is attributed to 3D induced effects and the specific airfoil pitching moment derivatives ($C_{m\text{foil}}$) present in the XFOIL polars.

## Verification
The proximity of both results to the 25-27% range confirms that the NLLT model and the MAC reference logic are consistent and provide a reliable baseline for the whole-aircraft stability analysis.
