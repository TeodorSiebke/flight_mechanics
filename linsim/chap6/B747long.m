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

g=9.8
W=2800000
m=W/g
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

As1=[Zw/(m-Zwdot) (Zq+m*u0)/(m-Zwdot)
(Mw+Mwdot*Zw/(m-Zwdot))/Iy (Mq+Mwdot*(Zq+m*u0)/(m-Zwdot))/Iy]

As2=[ Zw/(m) (Zq+m*u0)/(m)
(Mw+Mwdot*Zw/m)/Iy (Mq+Mwdot*(Zq+m*u0)/m)/Iy]

% Solve eigenvalue problem

[V,D]=eig(A)

% Non-dimensional eigenvectors x=( delta_u w q delta_tet)

vscal=[1/u0 1/u0 cbar/(2*u0) 1];
for j=1:4
v(:,j)=V(:,j).*vscal';
v(:,j)=v(:,j)/norm(v(:,j));
end

% State space model

%B=[0 0 0 0]';
%C=[1 0 0 0];
%D=[];
%sys=ss(A,B,C,D);
%initial(sys,[ 1 0 0 0],100)
B=[0 0 0 0]';
C=[0 0 1 0
   1 0 0 0];
D=[];
sys=ss(A,B,C,D);
initial(sys,[ 0 0 1 0],20)
%[y t x]=initial(sys,[ 0 0 1 0],200);
