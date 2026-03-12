% diag_30ms.m
addpath('flightsim');
fsm = make_fsim();
x_tr = [30 0.5 0 0.05 0 500 0 0.06 0.1]';
Kp = 0.01; Ki = 0.0001; Kd = 0.015; Kq = 0.4; h_target = 500;
xv = [x_tr(1:4); x_tr(6); 0]; % u, w, q, theta, h, h_int

try
    h_err = h_target - xv(5);
    h_dot = - ( -sin(xv(4))*xv(1) + cos(xv(4))*xv(2) );
    de = x_tr(8) - (Kp * h_err + Ki * xv(6) - Kd * h_dot) + Kq * xv(3);
    xf = [xv(1:4); 0; xv(5); 0; de; x_tr(9)];
    fprintf('xf size: %d x %d\n', size(xf,1), size(xf,2));
    xd = fplmod(0, xf, fsm);
    fprintf('xd size: %d x %d\n', size(xd,1), size(xd,2));
    res = [xd; h_err];
    fprintf('res size: %d x %d\n', size(res,1), size(res,2));
    fprintf('Success!\n');
catch ME
    fprintf('Error: %s\n', ME.message);
end
