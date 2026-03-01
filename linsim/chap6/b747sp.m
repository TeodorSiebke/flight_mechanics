% Simplified

%
% Boeing 747, p 165 Etkin
%
% State vector is x = (delta_u  w  q  delta_theta)
%

% US units -----------------------------------------------------
% Data in -units

% Aerodynamics

Xu=-2000
Xw=4000
Zu=-26000
Zw=-90000
Zq=-450000
Zwdot=2000
Mu=16000
Mw=-160000
Mq=-15000000
Mwdot=-17000

cbar=8.3

% Mass properties

g=10
m=280000
Iy=0.45*10^8

% Flight condition

u0=236
theta0=0

% The system matrix


A(1,1)=Xu/m
A(1,2)=Xw/m
A(1,3)=0
A(1,4)=-g*cos(theta0)
A(2,1)=Zu/(m-Zwdot)
A(2,2)=Zw/(m-Zwdot)
A(2,3)=(Zq+m*u0)/(m-Zwdot)
A(2,4)=-m*g*sin(theta0)/(m-Zwdot)
A(3,1)=(Mu+Mwdot*Zu/(m-Zwdot))/Iy
A(3,2)=(Mw+Mwdot*Zw/(m-Zwdot))/Iy
A(3,3)=(Mq+Mwdot*(Zq+m*u0)/(m-Zwdot))/Iy
A(3,4)=-Mwdot*m*g*sin(theta0)/(Iy*(m-Zwdot))
A(4,1)=0
A(4,2)=0
A(4,3)=1
A(4,4)=0

   eig(A)
clear
%------------------------------------

Zw=-90000
Zq=-450000
Zwdot=2000
Mw=-160000
Mq=-15000000
Mwdot=-17000

% Mass properties

g=10
m=280000
Iy=0.45*10^8

% Flight condition

u0=236
theta0=0

As1=[Zw/(m-Zwdot) (Zq+m*u0)/(m-Zwdot)
(Mw+Mwdot*Zw/(m-Zwdot))/Iy (Mq+Mwdot*(Zq+m*u0)/(m-Zwdot))/Iy]

As2=[ Zw/(m) (Zq+m*u0)/(m)
(Mw+Mwdot*Zw/m)/Iy (Mq+Mwdot*(Zq+m*u0)/m)/Iy]

As3=[ Zw/(m) (Zq+m*u0)/(m)
      Mw/Iy     Mq/Iy ]

% Solve eigenvalue problem

eig(As1)
eig(As2)
eig(As3)

% Non-dimensional eigenvectors
