% B747 data, p 165 Etkin
b=59.64
m=2.83176E6/9.81
Ix=0.247E8
Iy=0.449E8
Iz=0.673E8
Ixz=-0.212E7 % compare book
Ixz=1.315E6 % convert 9.7E5 slug*ft*ft

% US units
%Ix=1.82E7
%Iy=3.31E7
%Iz=4.97E7
%Ixz=9.70E5

II=[Ix 0 -Ixz
    0  Iy 0
    -Ixz 0 Iz]
%
   [V D]=eig(II)
%
% max (x^T II x^T)/(x^T x)
%
rg=sqrt(diag(D)/m)
