% Lateral B747 with control
% M=0.8 at 40000 ft altitude
u0=774;
theta0=0;
b=195.7;
% x = [ v p r phi]
A =[   -0.0558         0.0      -774.0000      32.2000
       -0.003865      -0.4342      0.4136       0.0
        0.001086      -0.006112   -0.1458       0.0
        0.0            1.0000      0.0          0.0    ];
% c = [da dr]'
B=[ 0.0       5.642
   -0.1431    0.1144
    0.003741 -0.4859
    0.0       0.0];
% psidot=r*sec(theta0)
% ydot=u0*psi*cos(theta0)+v

% Extend A to also include psi and y

Aext=[      A                  zeros(4,2)
   0 0 sec(theta0) 0         0               0
   1 0   0         0    u0*cos(theta0)       0];
Bext=[B;zeros(2,2)];
% Define output matrix with [v p r phi psi y]

Cext=[1 0 0 0 0 0
      0 1 0 0 0 0
      0 0 1 0 0 0
      0 0 0 1 0 0
      0 0 0 0 1 0
      0 0 0 0 0 1 ];
Dext=zeros(6,2);

vscal=[1/u0 b/(2*u0)  b/(2*u0) 1 1 ];

% Figures
close all
sys0=ss(Aext,(-15*pi/180)*Bext(:,1),Cext(1,:),Dext(1,1));
bode(sys0)
temp=input('<cr> for next');
clf

T=0:.01:30;
sys1=ss(Aext,(-15*pi/180)*Bext(:,1),Cext(1,:),Dext(1,1));
step(sys1,T);
title('Response to aileron -15 deg, Figure 7.30 (a)')
xlabel(' Time (s)');ylabel('v (fps)');
temp=input('<cr> for next');

T=0:.01:30;
sys2=ss(Aext,(-15*pi/180)*Bext(:,1),Cext(2:3,:),Dext(2:3,1));
step(sys2,T);
title('Response (p,r) to aileron -15 deg, Figure 7.30 (b)')
xlabel(' Time (s)');
temp=input('<cr> for next');

T=0:.01:30;
sys3=ss(Aext,(-15*pi/180)*Bext(:,1),(180/pi)*Cext(4:5,:),Dext(4:5,1));
step(sys3,T);
title('Response (\phi,\psi) to aileron -15 deg, Figure 7.30 (c)')
xlabel(' Time (s)');

% Use impulse instead...

