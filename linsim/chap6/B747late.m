% Lateral B747 with control
% M=0.8 at 40000 ft altitude
u0=774;
theta0=0;
b=195.7;
cbar=27.31;

% x = [ v p r phi]
A =[   -0.0558         0.0      -774.0000      32.2000
       -0.003865      -0.4342      0.4136       0.0
        0.001086      -0.006112   -0.1458       0.0
        0.0            1.0000      0.0          0.0    ];

vscal=[1/u0 b/(2*u0)  b/(2*u0) 1 ];

[V D]=eig(A)



for j=1:4
v(:,j)=V(:,j).*vscal';
v(:,j)=v(:,j)/norm(v(:,j));
end

% psidot=r*sec(theta0)
% ydot=u0*psi*cos(theta0)+v

% Extend A to also include psi and y

Aext=[      A                  zeros(4,2)
   0 0 sec(theta0) 0         0               0
   1 0   0         0    u0*cos(theta0)       0];

vscal=[1/u0 b/(2*u0)  b/(2*u0) 1 ];

vscalext=[vscal 1 1/b ];

[Vext Dext]=eig(Aext)

for j=1:6
vext(:,j)=Vext(:,j).*vscalext';
vext(:,j)=vext(:,j)/norm(vext(:,j));
end
