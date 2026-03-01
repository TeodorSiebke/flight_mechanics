%
% Run lifting-line code on the fuselage wind tunnel model wing
% This only the 600 mm span wing on the fuselage wind tunnel model
% This to be subtracted to get delta effect of fuselage

clear all;
close all;

% Assemble the model wing

acg = make_windexfusewing();

% Mass and center of gravity

acg.mass=0.0;
acg.cog=0.025;

% Define a single flight condition 

uoo = 35; 
AoA = 6*pi/180;
fs = flight_state(acg, uoo, 0.0, AoA, 0.0);

fs.delta_flap=0; % Dummy variable, implement Xfoil data before use

% Angles of attack

valpha = linspace(-5, 15, 21)' * pi/180; % '

% Compute an un-trimmed polar with 20 points

[vCL, vCD, vCM] = plainpolar(acg, fs, valpha);

% Approximate position of the aerodynamic center 

CLa = (vCL(10) - vCL(5)) / (valpha(10) - valpha(5))
CMa = (vCM(10) - vCM(5)) / (valpha(10) - valpha(5))
xac = acg.pref(1) - CMa/CLa * acg.cref
vCD

figure(1);
plot(valpha*180/pi, vCL, valpha*180/pi, vCM)
xlabel('Angle of attack [degree]');
ylabel('C_L and C_M');
grid on;

figure(2);
plot3(acg.pa(:,1),acg.pa(:,2),acg.pa(:,3),'kx', ...
      acg.pb(:,1),acg.pb(:,2),acg.pb(:,3),'ko')
  axis('equal')

