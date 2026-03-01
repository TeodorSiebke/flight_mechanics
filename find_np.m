% find_np.m
addpath('flightsim');
fsm = make_fsim();
vset = 30;
x = [vset 0.4 0.0 0.05 0 1000 0 0.0 0.1]'; 

% Trim
ivar = [1 2 4 8 9]; ifun = [1 2 3 6 8];
xtrim = x(ivar);
for iter = 1:50
    x(ivar) = xtrim;
    [xdot] = fplmod(0, x, fsm);
    xstep = 1e-7; J = zeros(5,5);
    for j = 1:5
        xh = x; xh(ivar(j)) = xh(ivar(j)) + xstep;
        xmh = x; xmh(ivar(j)) = xmh(ivar(j)) - xstep;
        [xdoth] = fplmod(0, xh, fsm);
        [xdotmh] = fplmod(0, xmh, fsm);
        J(:,j) = (xdoth(ifun) - xdotmh(ifun)) / (2*xstep);
    end
    ftrim = xdot(ifun); ftrim(5) = ftrim(5) - vset;
    if norm(ftrim) < 1e-9, break; end
    xtrim = xtrim - 0.7 * (J \ ftrim);
end
x(ivar) = xtrim;

% Derivatives w.r.t alpha (or w)
xstep = 1e-5;
xh = x; xh(2) = xh(2) + xstep;
xmh = x; xmh(2) = xmh(2) - xstep;
[xdoth] = fplmod(0, xh, fsm);
[xdotmh] = fplmod(0, xmh, fsm);

% dwdot/dw is roughly related to -Za or CL_alpha
% dqdot/dw is related to Ma or Cm_alpha
% qdot = (1/Iy) * (M_total)
% wdot = (1/m) * (Z_total) + q*u...
% Actually, let's just use Cm_alpha and CL_alpha directly if we can, 
% or compute them from forces.

% Easier: Compute forces and moments directly at two Alphas
[~, ~, info] = fplmod(0, x, fsm);
alpha1 = info.alpha;
Cm1 = info.Cm;
CL1 = info.CL;

x_h = x; x_h(2) = x_h(2) + 0.1; % Change w to change alpha
[~, ~, info_h] = fplmod(0, x_h, fsm);
alpha2 = info_h.alpha;
Cm2 = info_h.Cm;
CL2 = info_h.CL;

dCm_da = (Cm2 - Cm1) / (alpha2 - alpha1);
dCL_da = (CL2 - CL1) / (alpha2 - alpha1);

sm_nom = -dCm_da / dCL_da;
x_np = x(14) + sm_nom * fsm.cref; % Wait, x_cg is in fsm.cog(1)
x_cg_nom = fsm.cog(1);
x_np = x_cg_nom + sm_nom * fsm.cref;

fprintf('Nominal SM: %.2f%%\n', sm_nom * 100);
fprintf('Neutral Point x_np: %.4f m\n', x_np);
fprintf('cref: %.4f m\n', fsm.cref);

% Targets
sm_targets = [0.03, sm_nom, 0.15];
xcgs = x_np - sm_targets * fsm.cref;
fprintf('CG for 3%% SM:  %.4f m\n', xcgs(1));
fprintf('CG for Nom SM: %.4f m\n', xcgs(2));
fprintf('CG for 15%% SM: %.4f m\n', xcgs(3));
