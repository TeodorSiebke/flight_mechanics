% Boeing 747 Longitudinal example from Etkin

close all

theta0 = 0.0;
u0     = 774.0;

A = [ -0.00686854   0.01394951    0.0         -32.2  
      -0.09052721  -0.31506319  773.97653828    0.0
       0.00011865  -0.00102552   -0.42843608    0.0
       0.0          0.0           1.0           0.0  ];

B = [   -0.0002    9.6600
       -18.0800         0
        -1.1577         0
         0.0            0 ];

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
D =zeros(4,2);


%-----------------------------------------------------------------
% Elevator control input (first column of B)

sys=ss(A,B(:,1),(pi/180)*C(1,:),D(1,1));

bode(sys)
temp=input('<cr> for next');

% This sequence gives Figures 7.19 and 7.20
clf
T=0:.01:10;
sys=ss(A,B(:,1),(pi/180)*C(1,:),D(1,1));
step(sys,T);
title('Response to elevator (1 deg), Figure 7.19 (a)')
xlabel('Time (s)');ylabel('\Delta u (fps)');
temp=input('<cr> for next');

T=0:.01:60;
sys=ss(A,B(:,1),(pi/180)*C(2,:),D(2,1));
step(sys,T);
title('Response to elevator (1 deg), Figure 7.19 (b)')
xlabel('Time (s)');ylabel('\Delta \alpha (fps)');
temp=input('<cr> for next');

T=0:.01:10;
sys=ss(A,B(:,1),(pi/180)*C(3,:),D(3,1));
step(sys,T);
title('Response to elevator (1 deg), Figure 7.19 (c)')
xlabel('Time (s)');ylabel('\Delta \gamma (rad)');
temp=input('<cr> for next');

T=0:1:600;
sys=ss(A,B(:,1),(pi/180)*C(1,:),D(1,1));
step(sys,T);
title('Response to elevator (1 deg), Figure 7.20 (a)')
xlabel('Time (s)');ylabel('\Delta u (fps)');
temp=input('<cr> for next');

T=0:1:600;
sys=ss(A,B(:,1),(pi/180)*C(2,:),D(2,1));
step(sys,T);
title('Response to elevator (1 deg), Figure 7.20 (b)')
xlabel('Time (s)');ylabel('\Delta \alpha (fps)');
temp=input('<cr> for next');

T=0:1:600;
sys=ss(A,B(:,1),(pi/180)*C(3,:),D(3,1));
step(sys,T);
title('Response to elevator (1 deg), Figure 7.20 (c)')
xlabel('Time (s)');ylabel('\Delta \gamma (rad)');
temp=input('<cr> for next');

% Throttle control input (second column of B)

T=0:1:600;
sys=ss(A,B(:,2),(1/6)*C(1,:),D(1,2)) ; 
step(sys,T) ; 
title('Response to throttle, Figure 7.21 (a)')
xlabel('Time (s)');ylabel('\Delta u (fps)');
temp=input('<cr> for next');

sys=ss(A,B(:,2),(1/6)*C(2,:),D(2,2)) ; 
step(sys,T) ; 
title('Response to throttle, Figure 7.21 (b)')
xlabel('Time (s)');ylabel('\Delta \alpha (fps)');
temp=input('<cr> for next');

sys=ss(A,B(:,2),(1/6)*C(3,:),D(3,2)) ; 
step(sys,T) ; 
title('Response to throttle, Figure 7.21 (c)')
xlabel('Time (s)');ylabel('\Delta \gamma (rad)');


