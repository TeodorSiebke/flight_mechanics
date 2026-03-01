clear; close all;
fsm=make_fsim();

% Call model for trimmed point

x0=[
   29.9436    % u (m/s)	     
    1.8387    % w (m/s)	     
    0.0	      % q (rad/s)    
    0.0613    % theta (rad)  
    0.0	      % Distance (m) 
 1000.0	      % alt (m)	     
    0.0	      % fuel (kg)    
   -0.0250    % de (rad)     
    0.0906    % dp            
    ];

[xdot0]=fplmod(0,x0,fsm)

% Example uses illustrated below

% Linearization

x=x0;
nx = length(x);
nf = length(xdot0);
J = zeros(nf, nx);

xstep=10^(-5);
for j=1:nx
    xh=x; xmh=x;
    xh(j)=xh(j)+xstep;
    xmh(j)=xmh(j)-xstep;
    [xdoth]=fplmod(0,xh,fsm);
    [xdotmh]=fplmod(0,xmh,fsm);
    J(:,j)=(xdoth-xdotmh)/(2*xstep);
end

% Trim loop
% -------------------------------------------------------------------------
% ivar: The "Knobs" (independent variables the solver can change)
% [1]: u (Forward Vel), [2]: w (Vertical Vel), [4]: theta (Pitch Angle), 
% [8]: delta_e (Elevator), [9]: delta_p (Throttle)
ivar=[1 2 4 8 9];  

% ifun: The "Targets" (physics constraints we want to satisfy)
% [1]: u_dot=0, [2]: w_dot=0, [3]: q_dot=0 (No accelerations)
% [6]: alt_dot=0 (Level flight), [8]: Airspeed = Target
ifun=[1 2 3 6 8];  

xtrim=[30 0 0 0 0]'; % Initial guess
vset=30.0; % True airspeed is set

% Iterate until convergence

for iter=1:10
    x(ivar)=xtrim;
    [xdot]=fplmod(0,x,fsm);
    
    J = zeros(nf, nx);
    for j=1:nx
        xh=x; xmh=x;
        xh(j)=xh(j)+xstep;
        xmh(j)=xmh(j)-xstep;
        [xdoth]=fplmod(0,xh,fsm);
        [xdotmh]=fplmod(0,xmh,fsm);
        J(:,j)=(xdoth-xdotmh)/(2*xstep);
    end
    
    ftrim=xdot(ifun);
ftrim(5)=ftrim(5)-vset;
Jtrim=J(ifun,ivar);
xtrim=xtrim-Jtrim\ftrim;
fprintf('|f| %e\n',norm(ftrim));

end % End of Newton iteration

fprintf('alpha %e %e deg\n',xdot(9),xdot(9)*180/pi);
fprintf('airspeed %e m/s %e km/h\n',xdot(8),xdot(8)*3.6);

ev=eig(J(1:4,1:4))
T=2*pi./imag(ev)
T12=log(2)./real(ev)

% Evaluate 3D model

x3d=zeros(16,1);
x3d(1)=x(1);
x3d(3)=x(2);
x3d(8)=x(4);
x3d(12)=x(6);
x3d(13)=x(8);
x3d(16)=x(9);
[xdot]=fplmod3d(0,x3d,fsm);

x=x3d;
nx = length(x);
nf = length(xdot);
J = zeros(nf, nx);

xstep=10^(-5);
for j=1:nx
    xh=x; xmh=x;
    xh(j)=xh(j)+xstep;
    xmh(j)=xmh(j)-xstep;
    [xdoth]=fplmod3d(0,xh,fsm);
    [xdotmh]=fplmod3d(0,xmh,fsm);
    J(:,j)=(xdoth-xdotmh)/(2*xstep);
end

% Longitudinal Modal Analysis
% [1]: u (Vel), [3]: w (Heave), [5]: q (Pitch Rate), [8]: theta (Pitch Angle)
ivar=[1 3 5 8];  
ifun=[1 3 5 8];  

Along=J(ifun,ivar)
ev=eig(Along);
T=2*pi./imag(ev)
T12=log(2)./real(ev)

% Lateral-Directional Modal Analysis
% [2]: v (Sideslip), [4]: p (Roll Rate), [6]: r (Yaw Rate), [7]: phi (Bank Angle)
ivar=[2 4 6 7];  
ifun=[2 4 6 7];  

Alate=J(ifun,ivar)
ev=eig(Alate);
T=2*pi./imag(ev)
T12=log(2)./real(ev)

