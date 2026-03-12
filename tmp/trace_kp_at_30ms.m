% trace_kp_at_30ms.m
addpath('flightsim');
fsm = make_fsim();
v_target = 30;
h_target = 500;

% Fixed Damping
Kd = 0.015;
Kq = 0.4;
Ki = 0.0001;

% Trim
x_tr = [v_target 0.5 0 0.05 0 h_target 0 0.06 0.1]';
ivar = [1 2 4 8 9]; ifun = [1 2 3 6 8]; xt = x_tr(ivar);
for iter = 1:50
    x_tr(ivar) = xt; xd = fplmod(0, x_tr, fsm);
    res = xd(ifun); res(5) = res(5) - v_target;
    if norm(res) < 1e-10, break; end
    xs = 1e-7; J = zeros(5,5);
    for k = 1:5
        xh = x_tr; xh(ivar(k)) = xh(ivar(k)) + xs;
        xmh = x_tr; xmh(ivar(k)) = xmh(ivar(k)) - xs;
        xdh = fplmod(0, xh, fsm); xdm = fplmod(0, xmh, fsm);
        J(:,k) = (xdh(ifun) - xdm(ifun)) / (2*xs);
    end
    xt = xt - 0.7 * (J \ res);
end
x_tr(ivar) = xt;

kp_vec = logspace(-4, -1, 10);
for kp = kp_vec
    st_cl = [1 2 3 4 6]; n_aug = 6; A = zeros(n_aug, n_aug); xs = 1e-5;
    cl_sys = @(xv) cl_dynamics_debug(xv, x_tr, fsm, kp, Ki, Kd, Kq, h_target);
    for j = 1:n_aug
        xh = [x_tr(st_cl); 0]; xh(j) = xh(j) + xs;
        xmh = [x_tr(st_cl); 0]; xmh(j) = xmh(j) - xs;
        A(:,j) = (cl_sys(xh) - cl_sys(xmh)) / (2*xs);
    end
    eg = eig(A);
    fprintf('Kp = %.4f, Max Real = %e\n', kp, max(real(eg)));
end

function xdot_aug = cl_dynamics_debug(x_vec, x_tr, fsm, Kp, Ki, Kd, Kq, h_target)
    u = x_vec(1); w = x_vec(2); q = x_vec(3); theta = x_vec(4); h = x_vec(5); h_int = x_vec(6);
    h_err = h_target - h;
    h_dot = - ( -sin(theta)*u + cos(theta)*w );
    de = x_tr(8) - (Kp * h_err + Ki * h_int - Kd * h_dot) + Kq * q;
    xf = [u; w; q; theta; 0; h; 0; de; x_tr(9)];
    xd = fplmod(0, xf, fsm);
    xdot_aug = [xd(1); xd(2); xd(3); xd(4); xd(6); h_err];
end
