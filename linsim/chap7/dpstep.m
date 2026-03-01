% Boeing 747 Longitudinal example from Etkin

close all

theta0 = 0.0;
u0     = 774.0;

A = [ -0.00686854   0.01394951    0.0         -32.2  
      -0.09052721  -0.31506319  773.97653828    0.0
       0.00011865  -0.00102552   -0.42843608    0.0
       0.0          0.0           1.0           0.0  ];

B = [   9.66
           0.0
           0.0
           0.0];

ev=eig(A);

for i=1:length(ev)
  if imag(ev(i))~=0,
    Tper(i)=2*pi/imag(ev(i));
    nhalf(i)=(log(2)/(2*pi))*abs(imag(ev(i)))/abs(real(ev(i)));
    thalf(i)=log(2)/abs(real(ev(i)));
  end
end

% State vector is x = (delta_u  w  q  delta_theta)
% Output vector y = [ delta_u  alfa  gamma q ]

C =[1   0     0  0  
    0  1/u0   0  0  
    0 -1/u0   0  1  
    0   0     1  0 ];
D =zeros(4,1);

T=0:1:20;
sys=ss(A,B,(1/6)*C(1,:),D(1,1)) ; 
step(sys,T) ; 
title('Figure 7.21 (a)')
temp=input('<cr> for next');

sys=ss(A,B,(1/6)*C(2,:),D(2,1)) ; 
step(sys,T) ; 
title('Figure 7.21 (b)')
temp=input('<cr> for next');

sys=ss(A,B,(1/6)*C(3,:),D(3,1)) ; 
step(sys,T) ; 
title('Figure 7.21 (c)')


