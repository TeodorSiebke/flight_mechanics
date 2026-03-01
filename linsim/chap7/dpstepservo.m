% Boeing 747 Longitudinal example from Etkin

close all

taup=3.0;

theta0 = 0.0;
u0     = 774.0;

Ahat = [ -0.00686854   0.01394951    0.0         -32.2    9.66
         -0.09052721  -0.31506319  773.97653828    0.0    0.0
          0.00011865  -0.00102552   -0.42843608    0.0    0.0
          0.0          0.0           1.0           0.0    0.0
	  0.0          0.0           0.0           0.0   -1/taup ];

Bhat = [   0.0
           0.0
           0.0
           0.0
           1/taup ];

ev=eig(Ahat);

for i=1:length(ev)
  if imag(ev(i))~=0,
    T(i)=2*pi/imag(ev(i));
    nhalf(i)=(log(2)/(2*pi))*abs(imag(ev(i)))/abs(real(ev(i)));
    thalf(i)=log(2)/abs(real(ev(i)));
  end
end

% State vector is x = (delta_u  w  q  delta_theta delta_p)
% Output vector y = [ delta_u  alfa  gamma q ]

Chat =[1   0     0  0  0
       0  1/u0   0  0  0
       0 -1/u0   0  1  0
       0   0     1  0  0];
Dhat =zeros(5,1);

T=0:1:20;
sys=ss(Ahat,Bhat,(1/6)*Chat(1,:),Dhat(1,1)) ; 
step(sys,T) ; 
title('Figure 7.21 (a) with servo model')
temp=input('<cr> for next');

sys=ss(Ahat,Bhat,(1/6)*Chat(2,:),Dhat(2,1)) ; 
step(sys,T) ; 
title('Figure 7.21 (b) with servo model')
temp=input('<cr> for next');

sys=ss(Ahat,Bhat,(1/6)*Chat(3,:),Dhat(3,1)) ; 
step(sys,T) ; 
title('Figure 7.21 (c) with servo model')
