%
% Run lifting-line code on the Windex
%

clear all;
close all;

% Assemble the Windex model with tail surfaces

acg = make_windex();

% Mass and center of gravity

acg.mass=300;
acg.cog=0.25;

% Define a single flight condition 

uoo = 30; 
AoA = 6*pi/180;
fs = flight_state(acg, uoo, 0.0, AoA, 0.0);

fs.delta_elevator=0*pi/180;
fs.delta_rudder=0;
fs.delta_flap=0; % Dummy variable, implement Xfoil data before use
fs.delta_aileron=0;

% Solve for the vortex strength gamma (not flightpath angle) at fs

fobj = @(y) residual(acg, fs, y);
opt = optimset('Jacobian', 'on');
x = fsolve(fobj, zeros(size(acg.pa,1),1), opt);
fs.gamma = x;

% Compute and print coefficient data at fs

[CL, CD, CM] = postprocess(acg, fs)

% Display current L/D ratio 

CL / CD

% Angles of attack

valpha = linspace(-5, 20, 20)' * pi/180; % '

% Initialize matrices for CL CM and CD

mCL=[]; mCM=[]; mCD=[];

% Loop over several elevator settings

devec=[-10 -8 -6 -4 -2 0  2 0]*pi/180;

for ide=1:length(devec)

fs.delta_elevator=devec(ide);

% Compute an un-trimmed polar with 20 points

[vCL, vCD, vCM] = plainpolar(acg, fs, valpha);

% Approximate position of the aerodynamic center 

CLa = (vCL(10) - vCL(5)) / (valpha(10) - valpha(5))
CMa = (vCM(10) - vCM(5)) / (valpha(10) - valpha(5))
xac = acg.pref(1) - CMa/CLa * acg.cref

mCL=[mCL vCL];mCD=[mCD vCD];mCM=[mCM vCM];

end % end elevator loop

figure(1);
plot(valpha*180/pi, mCL, valpha*180/pi, mCM)
xlabel('Angle of attack [degree]');
ylabel('C_L and C_M');
grid on;

figure(2);
plot3(acg.pa(:,1),acg.pa(:,2),acg.pa(:,3),'kx', ...
      acg.pb(:,1),acg.pb(:,2),acg.pb(:,3),'ko')
  axis('equal')

