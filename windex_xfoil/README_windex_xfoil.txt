Generating polars using Xfoil

Viscous analysis is a NONLINEAR problem!

Convergence will often fail at higher alpha and discontinuities due to flap
deflections.

There may be multiple solutions for given Re and Alpha, as will the
real physics also have!

Make multiple tries from different starting alphas, use INIT to reset
when convergence reapeatly fails.

Start with low alpha and increase in small increments.

You can edit the polar files and adjust.

Make sure the rows have increasing values for alpha.

---------------------------

Basic airfoil data in file windex1200_originalprofil.dat

Typical xfoil session I use:

xfoil
load windex1200_originalprofil.dat
oper
visc
pacc
iter
aseq

Use a few different Reynolds number for better modeling and chose
a suitable range of alpha

Make sure to check every computed polar by reading using the matlab script

import_xfoil_polar.m

Plot CL and CM versus alpha and also CD vs CL look for strange
behaviour and discontinuities.

In particular, make sure each data point represents a converged solution.

I increase the number of iterations to 200.

For the tail airfoil, you need to modify the airfoil to represent
a geometry with deflected control surface.  One geometry for each
value of the control surface setting.

This can be done using Xfoil:

XFOIL   c>  load wxtail.dat
Enter airfoil name   s>  wxtail
 XFOIL   c>  gdes
.GDES   c>  f

Enter flap hinge x location   r>  0.7

  Top    surface:  y =  0.0355     y/t = 1.0
  Bottom surface:  y = -0.0355     y/t = 0.0

Enter flap hinge y location (or 999 to specify y/t)   r>  999

Enter flap hinge relative y/t location   r>  0.5

 Flap hinge: x,y =  0.70000 -0.00000

Enter flap deflection in degrees (+ down)   r>  10

 Top breaks: x,y =    0.70606  0.03440      0.70606  0.03440
 Bot breaks: x,y =    0.70305 -0.03493      0.70907 -0.03387
 Max thickness =     0.150345  at x =   0.347
 Max camber    =     0.036051  at x =   0.707

.GDES   c>  x

 Current airfoil nodes set from buffer airfoil nodes (  99 )

.GDES   c>  

 XFOIL   c>  save

Enter output filename   s>  wxtail+10.dat

 XFOIL   c>
 