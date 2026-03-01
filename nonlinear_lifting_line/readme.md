---
title: Nonlinear Lifting Line Library
author: David Eller
date: January 2018 
...

# Overview

This document gives is meant as a short guidance to the non-linear lifting line library for Octave/Matlab that has been developed for the course in Aircraft Aerodynamics at KTH. The first section is aimed as a tutorial of how to use the library, while the second part gives a short summary of the theoretical background and some of the limitation of this rather simple implementation.

# Using the library

The library consists of a collection of functions which enable to write programs (scripts)that perform various aerodynamic analyses of simplified aircraft configurations. Most of these functions operate on two types of objects (_structs_) that agglomerate the information needed to perform a set of aerodynamic analyses. 

## Construction 

Usually, a program will first construct an object which contains **geometry data** that does not change, such as the location of vortex elements, the location of the aerodynamic reference point and airfoil data.

~~~~ {.octave .numberLines}
acg = make_twinprop();
~~~~

Lifting surfaces expressed as vortex segments can be added by the utility 

~~~~ {.octave}
[acg, idx] = add_segments(acg, pts, chord, twist)
~~~~

which is passed a number of points *pts* along the leading edge of the wing (or tail surface) along with a vector of local chord and twist values (in radian). Lists of points are always expressed as a $[n \times 3]$ matrix with the $(x,y,z)$ components in the columns of the matrix. The function adds the corresponding segments to the object *acg* and returns the modified version along with a set of indices for the newly created vortex segments. These indices can be used to keep track of which vortex element belongs to which surface.

Airfoil data needs to be provided in the form of a lookup function which the solver (see below) will call in order to obtain two-dimensional aerodynamic coefficients for the airfoil segments associated with vortex segments. The lookup function can have any name but must have the form (arguments and return values)

~~~~ {.octave}
[cz, cdp, cm, cza] = my_airfoil_lookup(acg, fs, lalfa, re)
~~~~

where *lalfa* is a vector of effective, local two-dimensional angles of attack (in radian) seen by the airfoils, and *re* is a vector of Reynolds numbers computed with the local airfoil chord. The lookup function returns normal force coefficients *Cz*, profile drag (i.e. two-dimensional airfoil drag) coefficient *Cdp*, pitch moment coefficient *Cm* expressed around the 25% chord line, and normal force slope *Cza* (the derivative of the normal force coefficient with respect to the angle of attack). This lookup function is passed to the solver by attaching its handle to the geometry object:

~~~~ {.octave}
acg.foil = @my_airfoil_lookup;
~~~~

The airfoil data lookup function will normally need much more data, such as spline breakpoints or lookup tables. Such data should be attached to the geometry data object *acg* where it can be accessed when called by the solver. One of the examples (*twinprop*) demonstrates how to use airfoil data computed by *XFoil* (polars) in such a lookup function. 

The second type of object needed is the **flight state**. This object contains scalar data such as the velocity, angle of attack, sideslip angle, control surface deflections and rates of rotation, along with the effective velocity of the vortex segments. A flight state object also stores the distribution of vortex strengths that belongs to a particular flight condition.

~~~~ {.octave .numberLines startFrom="2"}
uoo = 250/3.6;  
alpha = 3.0*pi/180;
beta = 0.5*pi/180;
alt = 3000;
omega = [0 0 0];
fs = flight_state(acg, uoo, alt, alpha, beta, omega);

% set control variables 
fs.delta_flap = 0;
fs.delta_aileron = 7*pi/180;
fs.delta_elevator = 0;
fs.delta_rudder = 0;
~~~~

## Solving

Most library functions then need both an aircraft geometry object (*acg* above) and a flight state description object (*fs* above). An example is the simplest solver function 

~~~~ {.octave .numberLines startFrom="14"}
fs.gamma = pointsolve(acg, fs);
~~~~

which computes the vortex strength distribution **gamma** that corresponds to the flight state **fs** in the arguments. Once the object **fs** has been completed by *pointsolve()*, it can be used in the post-processing functions which compute forces, moments, and coefficients.

More advanced functions, such as trimming control surface deflections for a particular flight condition, can be expressed in a similar manner. In principle, there are two approaches to such higher-level problems:

1. Simultaneous trim;
2. Two-stage solution.

A simultaneous approach augments the vector of variables using in the non-linear solution procedure with the trim variables (such as control surface deflections) and the error vector with moment and force errors. The two-stage solution, on the other hand, solves for the vortex strengths $\Gamma$ repeatedly in an outer loop until the trim conditions are fulfilled. The first approach can be made more efficient, while the two-stage solution is much easier to program.

## Post-processing 

There are a number of functions that compute data derived from the vorticity distribution $\Gamma$ once that has been stored in a flight state object. On of them returns the non-dimensional force and moment coefficients:

~~~~ {.octave}
[CL, CD, CC, Cm, Cl, Cn] = coefficients(acg, fs)
~~~~

Note that these coefficients 

1. Follow Etkin's notation only when the default coordinate system is used for the geometry object *acg*;
2. Use the point *acg.pref* as reference for aerodynamic moments;
3. The drag coefficient *CD* includes the contribution of the fuselage drag.

Internally, this function calls a lower-level function 

~~~~ {.octave}
[Fsiv, Fsdp, Ms] = forces(acg, fs)
~~~~ 

which computes the force and moment contributions of each vortex segment. In the above call, *Fsiv* are the forces occurring due to inviscid flow (i.e. lift and induced drag) while *Fsdp* are the profile drag components. However, because the forces are returned as vectors with one entry for each vortex segment, any additional drag from fuselage or other non-lifting components are not included. Finally, *Ms* is a vector of the moment contributions. There is at the moment no model for the (usually small) aerodynamic moments generated by the fuselage.

## Example

The following example demonstrates how to solve for the vorticity $\Gamma$ corresponding to a certain flight condition and an aileron deflection of 7 degree. Using this distribution, the vortex forces and moments are determined and the spanwise load distribution plotted for the main wing in terms of local normal force coefficients. 

~~~~ {.octave .numberLines startFrom="15"}
% set control variables 
fs.delta_flap = 0;
fs.delta_aileron = 7*pi/180;
fs.delta_elevator = 0;
fs.delta_rudder = 0;

fs.gamma = pointsolve(acg, fs);
[Fsi, Fsv, Ms] = forces(acg, fs);
Fs = Fsi + Fsv;

% shortcuts for the index sets that tell us which 
% vortex belongs to which part 
irw = acg.vxwing;

% convert from vortex segment forces to local coefficients
Fzw = Fs(irw, 3);
y = 0.5*(acg.pa(irw,2) + acg.pb(irw,2));
As = abs(acg.pb(irw,2) - acg.pa(irw,2)) .* acg.chord(irw);
qoo = 0.5*fs.rho * fs.uref^2;
Czs = Fzw ./ (As * qoo);

plot(y, Czs);
xlabel('y [m]');
ylabel('c_z [-]');
title('Spanwise Load Distribution, d_a = 7deg');
~~~~

# Method

A lifting line method accounts only for lifting surfaces (no fuselage bodies or propulsion elements) and models each of the surfaces such as wings, canard or tail as a collection of vortex elements. Any number of vortex elements can be used to combine an arbitrary configuration of multiple lifting surfaces. The influence of a fuselage or other non-lifting bodies is not included. 

Each vortex consists of a straight bound element positioned at the 25% chord line, and two trailing vortex lines that extend downstream from the end points of the bound vortex towards infinity. In the default coordinate system, the downstream vortex direction is +x (this is a simplification for large angles of attack).  

A vortex strength $\Gamma$ is associated with each vortex segment; this strength determines the local circulation. A bound vortex with positive circulation generates a resulting force, the direction of which depends on the direction of the bound vortex line and the local velocity. Furthermore, the vortex generates an induced velocity field, that is, its presence changes the effective velocity everywhere to some degree. Far away from the vortex, this effect is weak and increases when approaching the vortex or its trailing segments. 

Hence, the set of all vortex elements and their strengths represented by a vector containing the values $\Gamma$ gives rise to a certain distribution of induced velocities $v_i(\boldsymbol{x})$ at some point $\boldsymbol{x}$. At each point, the motion of the body and the induced velocities results in an effective local velocity, so that a local angle of attack can be defined for each vortex segment. Additionally, the combination of strengths $\Gamma$ and local velocity results in a force $F_i$ acting on each vortex $i$. 

The core of the non-linear lifting-line method is to solve for some distribution of vortex strengths $\Gamma$ such that the resulting force generated by each vortex segment matches the force predicted by airfoil data for the local angle of attack given by body motion and the induced velocity at the center of the vortex segment. This is achieved by defining a residual, or error, vector $r_i(\Gamma) = F_i(\Gamma) - q_{\infty} S_i c_z(\alpha_i(\Gamma))$ for each vortex and solving for $\boldsymbol(r) = 0$.  

Being based on linearized potential flow, a linear lifting-line approach in itself is fairly limited. By including the effect of nonlinear airfoil characteristics, it becomes far more useful. Even though non-linear effects are only accounted for in the locally two-dimensional airfoil flow, such relevant quantities as control surface authority (and hence controllability and handling) and even an estimate of minimum speed and stall behavior can be obtained in principle, provided that accurate airfoil data is available.

Some limitations of the classical (linear) lifting-line method carry over to the non-linear version. Modelling a lifting surface as a set of vortex segments is only valid if there is no significant amount of span-wise flow. For this reason, only wings with sufficiently high aspect ratio (approximately $\Lambda > 5$) should be simplified in this way. Furthermore, flight conditions which incur considerable span-wise flow are unlikely to be represented accurately. Delta wings, for example, cannot be represented with any lifting-line approach in a meaningful manner. Similarly, large-scale separations are often (but not always!) accompanied by cross-flow; the stall behavior of a straight wing of large-aspect ratio may therefore be predicted reasonably well by a non-linear lifting-line method, while that is probably not the case for a significantly swept wing, where approaching stall leads to increasing amounts of span-wise flow in the boundary layer.

In principle, it is possible to approximately account for some effects of compressibility even with a simple lifting-line approach, by means of the Prandtl-Glauert transformation. That improvement is, however, not implemented in this library, which is hence limited to incompressible flow conditions (about $M < 0.3$).     


