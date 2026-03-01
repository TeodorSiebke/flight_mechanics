% verify_aero_change.m
% Quantifies the difference between 30m/s and 70m/s aero data in make_fsim.m

fsm = make_fsim();
alpha_idx = 6; % Alpha = 0 deg
de_idx = 2;    % Elevator = 0 deg

cl30 = fsm.cldat_30(alpha_idx, de_idx);
cl70 = fsm.cldat_70(alpha_idx, de_idx);

cd30 = fsm.cddat_30(alpha_idx, de_idx);
cd70 = fsm.cddat_70(alpha_idx, de_idx);

cm30 = fsm.cmdat_30(alpha_idx, de_idx);
cm70 = fsm.cmdat_70(alpha_idx, de_idx);

fprintf('--- Aero Data Comparison (Alpha=0, DE=0) ---\n');
fprintf('Parameter | 30 m/s  | 70 m/s  | Change %%\n');
fprintf('----------|---------|---------|----------\n');
fprintf('CL        | %7.4f | %7.4f | %+6.2f%%\n', cl30, cl70, (cl70-cl30)/abs(cl30)*100);
fprintf('CD        | %7.4f | %7.4f | %+6.2f%%\n', cd30, cd70, (cd70-cd30)/abs(cd30)*100);
fprintf('Cm        | %7.4f | %7.4f | %+6.2f%%\n', cm30, cm70, (cm70-cm30)/abs(cm30)*100);

% Test lon_aero interpolation
[cd_mid, cl_mid, cm_mid] = lon_aero(0, 0, 0, fsm, 50);

% Manual interpolation check for CL at Alpha=0, DE=0
cl_manual = 0.5 * fsm.cldat_30(6, 2) + 0.5 * fsm.cldat_50(6, 2);

fprintf('\n--- Interpolation Check (V=50 m/s) ---\n');
fprintf('Parameter | Expected | Actual  | Error\n');
fprintf('----------|----------|---------|----------\n');
fprintf('CL        | %8.4f | %8.4f | %8.4e\n', cl_manual, cl_mid, cl_mid - cl_manual);

% Check 30 m/s and 70 m/s results directly from lon_aero
[~, cl30_sim] = lon_aero(0, 0, 0, fsm, 30);
[~, cl70_sim] = lon_aero(0, 0, 0, fsm, 70);

fprintf('\n--- System Check (Direct lon_aero calls) ---\n');
fprintf('V       | Expected | Actual\n');
fprintf('--------|----------|----------\n');
fprintf('30 m/s  | %8.4f | %8.4f\n', fsm.cldat_30(6,2), cl30_sim);
fprintf('70 m/s  | %8.4f | %8.4f\n', fsm.cldat_70(6,2), cl70_sim);
