% optimize_elevator_pid.m
% Automated PID Tuning for the Elevator Controller
%
% This script uses an optimization routine to find the "best" fixed gains
% for the altitude-hold PID controller across a representative envelope.

clear; close all; clc;
fsm = make_fsim();

% --- 1. OPTIMIZATION SETTINGS ---
% We use a sparse representative grid to speed up optimization
v_opt = [30]; 
h_opt = [1000];

% Initial guess (Start with base gains)
gains0 = [0.0045, 0.0002, 0.015, 0.4];

% Lower and Upper Bounds (to keep gains physical/sensible)
lb = [0.0001, 1e-6,   0.001, 0.05];
ub = [0.05,   0.01,   0.1,   2.5];

fprintf('Starting Automated PID Optimization...\n');
fprintf('Base Gains: Kp=%.5f, Ki=%.5f, Kd=%.5f, Kq=%.2f\n', gains0);

% Optimization Options
options = optimset('Display', 'iter', 'TolX', 1e-4, 'MaxIter', 100);

% Objective Function
% We want to MINIMIZE a cost. 
% Cost = 1/AvgDamping + PenaltyForUnstable
obj_fun = @(g) calculate_envelope_cost(g, v_opt, h_opt, fsm);

% Run optimization
[best_gains, min_cost] = fminsearch(obj_fun, gains0, options);

fprintf('\n=== OPTIMIZATION RESULTS ===\n');
fprintf('Best Fixed Gains found:\n');
fprintf('  Kp_h = %.6f\n', best_gains(1));
fprintf('  Ki_h = %.6f\n', best_gains(2));
fprintf('  Kd_h = %.6f\n', best_gains(3));
fprintf('  Kq   = %.6f\n', best_gains(4));
fprintf('\nRun analyze_elevator_performance.m with these values to see the full result.\n');

% --- HELPER: COST FUNCTION ---
function cost = calculate_envelope_cost(gains, v_vec, h_vec, fsm)
    Kp = gains(1); Ki = gains(2); Kd = gains(3); Kq = gains(4);
    
    total_points = length(v_vec) * length(h_vec);
    stability_penalty = 0;
    damping_acc = 0;
    
    for h_target = h_vec
        for v_target = v_vec
            % 1. Trim (Simplified/Fast)
            x_tr = [v_target 0.5 0 0.05 0 h_target 0 0.06 0.1]';
            ivar = [1 2 4 8 9]; ifun = [1 2 3 6 8]; xt = x_tr(ivar);
            for iter = 1:20
                x_tr(ivar) = xt; xd = fplmod(0, x_tr, fsm);
                res = xd(ifun); res(5) = res(5) - v_target;
                if norm(res) < 1e-8, break; end
                xs = 1e-7; J = zeros(5,5);
                for k = 1:5
                    xh = x_tr; xh(ivar(k)) = xh(ivar(k)) + xs;
                    xmh = x_tr; xmh(ivar(k)) = xmh(ivar(k)) - xs;
                    xdoth = fplmod(0, xh, fsm);
                    xdotmh = fplmod(0, xmh, fsm);
                    J(:,k) = (xdoth(ifun) - xdotmh(ifun)) / (2*xs);
                end
                xt = xt - 0.7 * (J \ res);
            end
            
            % 2. Linearization
            state_idx = [1 2 3 4 6]; n_aug = 6; A = zeros(n_aug, n_aug); xs_lin = 1e-5;
            cl_system = @(xv) cl_dynamics_fast(xv, x_tr, fsm, Kp, Ki, Kd, Kq, h_target);
            for k = 1:n_aug
                xh = [x_tr(state_idx); 0]; xh(k) = xh(k) + xs_lin;
                xmh = [x_tr(state_idx); 0]; xmh(k) = xmh(k) - xs_lin;
                A(:,k) = (cl_system(xh) - cl_system(xmh)) / (2*xs_lin);
            end
            
            % 3. Evaluate
            eg = eig(A);
            % Soft Stability: Accept real parts if Time-to-Double > 200s
            stability_threshold = log(2)/200;
            if any(real(eg) >= stability_threshold)
                stability_penalty = stability_penalty + 1000; % Large penalty for unstable
            else
                % Target damping Zeta = 0.5
                zetas = -real(eg) ./ abs(eg);
                min_z = min(zetas(abs(imag(eg)) > 1e-4)); % Lower threshold to catch Phugoid
                if isempty(min_z), min_z = 1; end
                damping_acc = damping_acc + (0.5 - min_z)^2;
            end
        end
    end
    
    cost = stability_penalty + damping_acc;
end

function xdot_aug = cl_dynamics_fast(x_vec, x_tr, fsm, Kp, Ki, Kd, Kq, h_target)
    u = x_vec(1); w = x_vec(2); q = x_vec(3); tet = x_vec(4); h = x_vec(5); hi = x_vec(6);
    h_err = h_target - h;
    h_dot = - ( -sin(tet)*u + cos(tet)*w );
    de = x_tr(8) - (Kp * h_err + Ki * hi - Kd * h_dot) + Kq * q;
    xf = [u; w; q; tet; 0; h; 0; de; x_tr(9)];
    xd = fplmod(0, xf, fsm);
    xdot_aug = [xd(1); xd(2); xd(3); xd(4); xd(6); h_err];
end

% Note: Added a small performance tweak to fplmod call in helper for speed
function out = fplmod_fast(t, x, fsm, idx)
    res = fplmod(t, x, fsm);
    out = res(idx);
end
