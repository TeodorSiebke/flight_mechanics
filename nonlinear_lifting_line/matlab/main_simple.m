%
% run lifting-line code on a simple example 
%

clear all;

% construct a very simple aircraft model with twisted wing and simple HTP
nw = 24;
acg = make_example(nw);

% For the second windtunnnel lab this aircraft model can be used instead, 
% but needs to be filled in with relevant geometric dimensions!
%acg = make_wtm2(nw);

% define a single flight condition 
uoo = 80; 
AoA = 6*pi/180;
fs = flight_state(acg, uoo, 0.0, AoA, 0.0);

% g0 = zeros(size(acg.pa,1),1);
% [r, J] = residual(acg, fs, g0);

% solve for the vortex strength gamma at fs
fobj = @(y) residual(acg, fs, y);
opt = optimset('Jacobian', 'on');
x = fsolve(fobj, zeros(size(acg.pa,1),1), opt);
fs.gamma = x;

% compute and print coefficient data at fs
[CL, CD, CM] = postprocess(acg, fs)

% display current L/D ratio 
CL / CD

% plot spanwise distributions
[cxyz, vt] = segment_coefficients(acg, fs, fs.gamma);
pmid = 0.5*(acg.pa + acg.pb);
yw = pmid(1:nw, 2);
yt = pmid((nw+1):end, 2);

% local normal force coefficient Cz on main wing and HTP
figure(1);
plot(yw, cxyz(1:nw,3), yt, cxyz((nw+1):end, 3));
xlabel('Span coordinate [m]');
ylabel('Section lift coefficient [-]');
legend('Wing', 'HTP');

% compute an un-trimmed polar with 20 points

valpha = linspace(-5, 20, 20)' * pi/180;
[vCL, vCD, vCM] = plainpolar(acg, fs, valpha);

% approximate position of the aerodynamic center 
CLa = (vCL(10) - vCL(1)) / (valpha(10) - valpha(1))
CMa = (vCM(10) - vCM(1)) / (valpha(10) - valpha(1))
xac = acg.pref(1) - CMa/CLa * acg.cref

% because the simple example aircraft uses an airfoil function which 
% just returns Cz = 2*pi*alpha, the simple example aircraft won't stall...
figure(2);
plot(valpha*180/pi, vCL, valpha*180/pi, 10*vCM);
xlabel('Angle of attack [degree]');
legend('C_L', '10*C_M');
grid on;

% again, for the simple model, this is very optimistic because there is no 
% drag contribution from flow separation at high AoA in this model - 
% improve by plugging in a better airfoil function, like the one used in the 
% model for the second windtunnel lab!

figure(3);
plot(vCL, vCL ./ vCD);
xlabel('Lift coefficient [-]');
ylabel('L/D');
