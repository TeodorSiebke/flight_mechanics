% Main program

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
    0.0    % dp (Throttle)          
    ];

   [xdot0]=fplmod(0,x0,fsm)

% Example uses illustrated below

% Linearization

x=x0;

% Trim loop

ivar=[1 2 4 8];  % Selects variables to change during trim iteration
ifun=[1 2 3 8];  % Selects functions that should be considered in trim

xtrim=[30 0 0 0]'; %' Initial guess

mfpl=[270 370]
vstart=[24 28]
for imas=1:2

fsm.mass=mfpl(imas)
vset=vstart(imas); % True airspeed is set

for isp=1:400

% Iterate until convergence

for iter=1:10

x(ivar)=xtrim;
[xdot]=fplmod(0,x,fsm);

xstep=10^(-5);
for j=1:length(x)
xh=x;
xmh=x;
xh(j)=xh(j)+xstep;
xmh(j)=xmh(j)-xstep;
[xdoth]=fplmod(0,xh,fsm);
[xdotmh]=fplmod(0,xmh,fsm);
J(:,j)=(xdoth-xdotmh)/(2*xstep);
end

ftrim=xdot(ifun);
ftrim(4)=ftrim(4)-vset;
Jtrim=J(ifun,ivar);
xtrim=xtrim-Jtrim\ftrim;
fprintf('|f| %e\n',norm(ftrim));

end % End of Newton iteration

fprintf('alpha %e %e deg\n',xdot(9),xdot(9)*180/pi);
fprintf('airspeed %e m/s %e km/h\n',xdot(8),xdot(8)*3.6);

ev=eig(J(1:4,1:4))
T=2*pi./imag(ev)
T12=log(2)./real(ev)

vvec(isp,imas)=vset;
sinkrate(isp,imas)=xdot(6);
devec(isp)=x(8);

vset=vset+.1
end % isp loop
end % imas

close all
plot(vvec(:,1)*3.6,sinkrate(:,1),vvec(:,2)*3.6,sinkrate(:,2))
legend('m = 270 kg','m = 370 kg')
xlabel('Air speed (km/h)')
ylabel('Sink rate (m/s)')
