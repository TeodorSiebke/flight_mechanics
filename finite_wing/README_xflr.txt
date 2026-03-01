liftline1.m CLa=4.1962
>> rectwing_polar1
>> plot(poldat(:,1),poldat(:,2),[0 10],[0 CL*pi/18])


xflr5

new project

File/xfoil direct analysis

Design/ naca foils

chose "12" for NACA0012

Batch analysis

Type 1
Range of Re
Analysis range (Alpha) -1 to 20 

b =    1.1250
chord =    0.2472
nu=1.5e-5
Relow=chord*10/nu
ans =   1.6481e+05
Rehigh=chord*50/nu
8.2407e+05

chose 100 000 to 900 000 step 100 000

File/wing and plane design/Wing-Plane

Define a New Wing

c=247.2
b/2=562.5

Analysis/Define an analysis

Free stream speed
Type 1
LLT

Polars
Current polar
export

----------

Convergence problem may depend on:

Not enough xfoil data covering the entire range of Reynolds number

RE_min= minimum speed * minimum chord / viscosity

RE_max= maximum speed * maximum chord / viscosity

You can also try to reduce the number of horse shoe vortices used.
Fewer elements means a less accurate solution, but the nonlinear
problem that needs to be solved is smaller and therefore usually
easier to solve.
