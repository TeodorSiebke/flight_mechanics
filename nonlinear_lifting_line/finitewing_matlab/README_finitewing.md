# Finite Wing NLLLT Analysis: Documentation

This folder contains the implementation of the **Non-Linear Lifting Line Theory (NLLLT)** for the finite wing model, adapted from the KTH Aircraft Aerodynamics library.

## Project Overview
The goal of this code is to compute the aerodynamic performance (Lift, Drag, Pitching Moment) of a finite wing using 2D airfoil data (from XFOIL) and compare it against experimental wind tunnel data.

## File Structure

| File | Description |
| :--- | :--- |
| `main_finitewing.m` | **Primary Entry Point**. Runs the NLLLT solver and generates comparison plots. |
| `make_finitewing.m` | Defines the wing geometry, solver mesh (courseness), and moment reference point. |
| `build_finitewing.m` | Dynamically loads XFOIL polar files and constructs the aerodynamic grid. |
| `finitewing_foils.m` | Wrapper function used by the solver to look up airfoil coefficients. |
| `check_xfoil_smoothing.m` | **Diagnostic Tool**. Run this to visualize your XFOIL polars and check for noise. |

## How to Run

1.  **Generate Polars**: Follow the XFOIL walkthrough to generate files named `h15_re*e5.txt`. 
2.  **Verify Data**: Run `check_xfoil_smoothing.m` in MATLAB to ensure the curves are smooth.
3.  **Run Analysis**: Run `main_finitewing.m`. 

## Configuration Guide

### Geometry and Reference Point
Open `make_finitewing.m` to adjust:
- **Courseness (`nw`)**: Number of spanwise segments (default=30, increase for better precision).
- **Reference Point (`acg.pref`)**: The $[x,y,z]$ location used for moment calculations (should match your wind tunnel balance point).
- **Planform**: The $p_{tip}$ and $p_{root}$ coordinates define the leading edge shape.

### Solver Controls
Open `main_finitewing.m` to adjust:
- **Alpha Limits**: `alpha_limit_min` and `alpha_limit_max` define the simulation range.
- **Pitch Offset**: `alpha_offset` (default -2.5) aligns the wind tunnel data with the XFOIL coordinate system.
- **Data Source**: Change `wt_file` to point to different wind tunnel runs (e.g., `u20.txt`, `u40.txt`).

## Technical Implementation Details
- **Continuation Method**: The solver uses the solution from the previous angle of attack as the initial guess for the next. This significantly improves stability and speed.
- **Convergence Guard**: The code strictly validates the `fsolve` exit flag. If the solver fails to converge at a specific angle, that point is discarded (set to `NaN`) to prevent misleading results in the plots.
- **Viscous Blending**: At extreme angles, the library blends XFOIL data into a flat-plate model to maintain numerical stability.
