% find_np_proper.m
addpath('flightsim');
fsm = make_fsim();

% Get reference values
alpha_nom = 0.05; % Approximate trim alpha
deltae = 0.0;
qbar = 0.0;

[CD1, CL1, Cm1] = lon_aero(alpha_nom, deltae, qbar, fsm);
[CD2, CL2, Cm2] = lon_aero(alpha_nom + 0.01, deltae, qbar, fsm);

dCL_da = (CL2 - CL1) / 0.01;
dCm_da = (Cm2 - Cm1) / 0.01;

% The Cm in lon_aero is about the pref (0.25m).
% fsm.cog(1) is 0.369m.
% Cm_cg = Cm_pref + CL * (x_cg - x_pref)/cref
% dCm_cg/da = dCm_pref/da + dCL/da * (x_cg - x_pref)/cref
% At Neutral Point, dCm_cg/da = 0
% dCm_pref/da + dCL/da * (x_np - x_pref)/cref = 0
% (x_np - x_pref)/cref = - (dCm_pref/da) / (dCL_da)
% x_np = x_pref - cref * (dCm_da/dCL_da)

x_ref = fsm.pref(1); % 0.25
c_ref = fsm.cref;     % 0.641

x_np = x_ref - c_ref * (dCm_da / dCL_da);

fprintf('Neutral Point x_np: %.4f m\n', x_np);

% Nominal SM
x_cg_nom = fsm.cog(1); % 0.369
sm_nom = (x_np - x_cg_nom) / c_ref;
fprintf('Nominal SM: %.2f%%\n', sm_nom * 100);

% Targets
sm_targets = [0.03, sm_nom, 0.15];
xcgs = x_np - sm_targets * c_ref;
fprintf('CG for 3%% SM:  %.4f m (SM=3%%)\n', xcgs(1));
fprintf('CG for Nom SM: %.4f m (SM=%.2f%%)\n', xcgs(2), sm_nom*100);
fprintf('CG for 15%% SM: %.4f m (SM=15%%)\n', xcgs(3));
