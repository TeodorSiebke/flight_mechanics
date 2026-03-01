%
% Very simple lifting line model based on Appendix E
% Approximate trapezoidal wing with simple rectangular planform
%
b=1.125;     % Total span (m)
bw=0.125;    % Width of constant chord mid section (m)
croot=0.277; % Root section chord (m)
ctip=0.21;   % Tip chord (m)
xstn=[.3577 0.3006 0.3006 .3577]; % Leading edge coordinates

% Approximate geometry (DO NOT USE)

span=b;                              % Actual span
Sref=(croot+ctip)*(b-bw)/2+croot*bw; % Actual projected wing area
chord=Sref/span;                     % Mean chord

N=30; % Number of segments

ctetm=chord; % Chord constant, modify to give real wing shape
alfa=1;      % Unit angle of attack
cla=2*pi;    % Section lift slope

for im=1:N
for in=1:N

vtetm=pi*im/(N+1);
A(im,in)=sin(in*vtetm)+ctetm*cla*in*sin(in*vtetm)/(4*span*sin(vtetm));
end
rhs(im,1)=ctetm*cla*alfa/(4*span);
end
x=A\rhs;
CL=x(1)*pi*span/chord
