function diag_unstable_mode()
    % diag_unstable_mode.m
    % Detailed eigenvalue analysis for Speed-Hold PID
    
    addpath('../flightsim');
    fsm = make_fsim();
    
    % Test Point: 25 m/s at 0m (Sea Level)
    v_target = 25;
    h_target = 0;
    
    % --- DIAGNOSTIC: 10000m @ 25 m/s ---
    v_target = 25; h_target = 10000;
    Kp = 0.001785; Ki = 0.000259; Kd = 0.018147; Kq = 0.054794;
    
    % 1. Trim
    x_tr = [v_target 0.5 0 0.05 0 h_target 0 0.06 0.1]';
    ivar = [1 2 4 8 9]; % u, w, theta, de, dp
    ifun = [1 2 3 6 8]; % udot, wdot, qdot, hdot, v-v_target
    xt = x_tr(ivar);
    for iter = 1:50
        x_tr(ivar) = xt;
        xd = fplmod(0, x_tr, fsm);
        res = xd(ifun); res(5) = res(1) - 0; % res(5) should be udot or vdot
        % Actually, standard trim: udot, wdot, qdot, hdot, v=v_target
        % But the logic in elevator_level_flight used res(5) = res(5) - v_target
        % Let's stick to a reliable trim:
        xd = fplmod(0, x_tr, fsm);
        res = [xd(1); xd(2); xd(3); xd(6); sqrt(x_tr(1)^2+x_tr(2)^2)-v_target];
        
        if norm(res) < 1e-10, break; end
        xs = 1e-8; J = zeros(5,5);
        for j = 1:5
            xh = x_tr; xh(ivar(j)) = xh(ivar(j)) + xs;
            xmh = x_tr; xmh(ivar(j)) = xmh(ivar(j)) - xs;
            xdh = fplmod(0, xh, fsm);
            xdm = fplmod(0, xmh, fsm);
            resh = [xdh(1); xdh(2); xdh(3); xdh(6); sqrt(xh(1)^2+xh(2)^2)-v_target];
            resm = [xdm(1); xdm(2); xdm(3); xdm(6); sqrt(xmh(1)^2+xmh(2)^2)-v_target];
            J(:,j) = (resh - resm) / (2*xs);
        end
        xt = xt - 0.7*(J \ res);
    end
    x_tr(ivar) = xt;
    
    % 2. Linearization
    % States: [u w q theta v_int]
    s_idx = [1 2 3 4]; 
    n = 5;
    A = zeros(n,n);
    xsl = 1e-5;
    
    for k = 1:n
        xh = [x_tr(s_idx); 0]; xh(k) = xh(k) + xsl;
        xmh = [x_tr(s_idx); 0]; xmh(k) = xmh(k) - xsl;
        
        xdh = speed_hold_cl_dynamics(xh, x_tr, fsm, Kp, Ki, Kd, Kq, v_target);
        xdm = speed_hold_cl_dynamics(xmh, x_tr, fsm, Kp, Ki, Kd, Kq, v_target);
        A(:,k) = (xdh - xdm) / (2*xsl);
    end
    
    [V, D] = eig(A);
    eg = diag(D);
    
    disp(['--- Stability Diagnostics for Speed-Hold (V=', num2str(v_target), ') ---']);
    states = {'u', 'w', 'q', 'theta', 'v_int'};
    
    for i = 1:length(eg)
        participation = abs(V(:,i)) / sum(abs(V(:,i))) * 100;
        [val, idx] = sort(participation, 'descend');
        p_str = sprintf('%s (%.1f%%), %s (%.1f%%)', states{idx(1)}, val(1), states{idx(2)}, val(2));
        
        fprintf('Mode %d: Eigenvalue = %.5f%+.5fi\n', i, real(eg(i)), imag(eg(i)));
        if real(eg(i)) > 0
            fprintf('  *** WARNING: UNSTABLE! (T_double = %.1fs) ***\n', log(2)/real(eg(i)));
        else
            fprintf('  Stable (T_half = %.1fs)\n', log(2)/abs(real(eg(i))));
        end
        fprintf('  Participation: %s\n\n', p_str);
    end
end

function xd_aug = speed_hold_cl_dynamics(x_vec, x_tr, fsm, Kp, Ki, Kd, Kq, v_target)
    u = x_vec(1); w = x_vec(2); q = x_vec(3); tet = x_vec(4); vi = x_vec(5);
    v_curr = sqrt(u^2 + w^2); v_err = v_target - v_curr;
    
    % Physics xdot for v_dot
    xf_tmp = [u; w; q; tet; 0; x_tr(6); 0; x_tr(8); x_tr(9)];
    xd_tmp = fplmod(0, xf_tmp, fsm);
    v_dot = (u*xd_tmp(1) + w*xd_tmp(2)) / v_curr;
    
    de = x_tr(8) + (Kp * v_err + Ki * vi - Kd * v_dot) + Kq * q;
    xf = [u; w; q; tet; 0; x_tr(6); 0; de; x_tr(9)];
    xd = fplmod(0, xf, fsm);
    xd_aug = [xd(1); xd(2); xd(3); xd(4); v_err];
end
