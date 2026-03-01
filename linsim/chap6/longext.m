load B747long

A=[ALONG zeros(4,2) 
  cos(theta0) sin(theta0) 0 -u0*sin(theta0)  0  0 
 -sin(theta0) cos(theta0) 0 -u0*cos(theta0)  0  0  ]

C=[0 0 0 0 1 0
   0 0 0 0 0 -1]

  sys=ss(A,zeros(6,1),C,zeros(2,1)) 

[y t x]=  initial(sys,[0 0.01 0 0 0 0],400);

plot(y(:,1),y(:,2))

%B=[]
  [V D]=eig(A);
cbar=27.31;
vscal=[1/u0 1/u0 cbar/(2*u0) 1 1 1 ];
for j=1:6
v(:,j)=V(:,j).*vscal';
v(:,j)=v(:,j)/norm(v(:,j));
end
