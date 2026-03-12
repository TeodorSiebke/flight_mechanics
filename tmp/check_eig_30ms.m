% check_eig_30ms.m
addpath('flightsim');
fsm = make_fsim();
v_target = 30;
h_target = 500;

% Gains
Kp = 0.01; Ki = 0.0001; Kd = 0.015; Kq = 0.4;

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

% Linearize
st_cl = [1 2 3 4 6]; n_aug = 6; A = zeros(n_aug, n_aug); xs = 1e-5;
cl_sys = @(xv) [fplmod(0, [xv(1:4); 0; xv(5); 0; x_tr(8)-(Kp*(h_target-xv(5))+Ki*xv(6)-Kd*(-( -sin(xv(4))*xv(1) + cos(xv(4))*xv(2)))); x_tr(9)], fsm); h_target-xv(5)];

for j = 1:6
    xh = [x_tr(1:4); x_tr(6); 0]; xh(j) = xh(j) + xs;
    xmh = [x_tr(1:4); x_tr(6); 0]; xmh(j) = xmh(j) - xs;
    resh = cl_sys(xh); resmh = cl_sys(xmh);
    A(:,j) = (resh([1 2 3 4 6 14]) - resmh([1 2 3 4 6 14])) / (2*xs);
end

fprintf('--- 30 m/s at 500m ---\n');
fprintf('Gains: Kp=%.4f, Ki=%.6f, Kd=%.4f, Kq=%.4f\n', Kp, Ki, Kd, Kq);
fprintf('Eigenvalues:\n');
eg = eig(A);
disp(eg);
fprintf('Max Real Part: %e\n', max(real(eg)));
