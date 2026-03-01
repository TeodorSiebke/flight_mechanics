load S2000-4.mat;Init2000;ev=eig(ALONG)
T=2*pi./imag(ev)
Thalf=log(2)./real(ev)

load S2000-1.mat;Init2000;ev=eig(ALATE);disp(ev);T=2*pi./imag(ev);disp(T);Thalf=log(2)./real(ev)

[V D]=eig(ALATE)

scalevec

V.*scalelate

ALONGext=[ALONG zeros(4,2) 
  cos(theta0) sin(theta0) 0 -u0*sin(theta0)  0  0 
 -sin(theta0) cos(theta0) 0 -u0*cos(theta0)  0  0  ]

ALATEext=[      ALATE            zeros(4,2)
   0 0 sec(theta0) 0         0               0
   1 0   0         0    u0*cos(theta0)       0];

scalelong=[1/u0 1/u0 CREF(1)/(2*u0) 1 1 1]'
scalelate=[1/u0 BREF(1)/(2*u0) BREF(1)/(2*u0) 1 1/BREF(1)]'
