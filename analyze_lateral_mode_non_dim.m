% analyze_lateral_mode_non_dim.m
% Non-dimensionalized eigenvector analysis for Lateral-Directional modes
addpath('flightsim');
fsm = make_fsim();
vset = 30; % Analyze at moderate speed
bref = fsm.bref;

% 1. Longitudinal Trim
x_lon = [vset 0.4 0.0 0.05 0 1000 0 0.0 0.1]'; 
ivar_lon = [1 2 4 8 9]; ifun_lon = [1 2 3 6 8];
xtrim = x_lon(ivar_lon);
for iter = 1:50 
    x_lon(ivar_lon) = xtrim;
    [xdot_lon] = fplmod(0, x_lon, fsm);
    xstep = 1e-7; J = zeros(5, 5);
    for j = 1:5
        xh = x_lon; xh(ivar_lon(j)) = xh(ivar_lon(j)) + xstep;
        xmh = x_lon; xmh(ivar_lon(j)) = xmh(ivar_lon(j)) - xstep;
        [xdoth] = fplmod(0, xh, fsm); [xdotmh] = fplmod(0, xmh, fsm);
        J(:,j) = (xdoth(ifun_lon) - xdotmh(ifun_lon)) / (2*xstep);
    end
    ftrim = xdot_lon(ifun_lon); ftrim(5) = ftrim(5) - vset;
    if norm(ftrim) < 1e-9, break; end
    xtrim = xtrim - 0.7 * (J \ ftrim);
end
x_lon(ivar_lon) = xtrim;

% 2. 3D Linearization
x3d = zeros(16,1);
x3d(1) = x_lon(1); x3d(3) = x_lon(2); x3d(5) = x_lon(3);
x3d(8) = x_lon(4); x3d(12) = x_lon(6); x3d(13) = x_lon(8); x3d(16) = x_lon(9);

lat_idx = [2 4 6 7]; % v, p, r, phi
Alat = zeros(4,4); xstep = 1e-5;
for j = 1:4
    idx = lat_idx(j);
    xh = x3d; xh(idx) = xh(idx) + xstep;
    xmh = x3d; xmh(idx) = xmh(idx) - xstep;
    [xdoth] = fplmod3d(0, xh, fsm); [xdotmh] = fplmod3d(0, xmh, fsm);
    Alat(:,j) = (xdoth(lat_idx) - xdotmh(lat_idx)) / (2*xstep);
end

% 3. Eigenvector Analysis
[V, D] = eig(Alat);
eg = diag(D);

% Non-dimensionalization scaling factors
% State: [v (m/s), p (rad/s), r (rad/s), phi (rad)]
% Scaled: [v/V, p*b/(2V), r*b/(2V), phi]
scaling = [1/vset, bref/(2*vset), bref/(2*vset), 1];

fprintf('\nLateral-Directional Mode Coupling Analysis at %.0f m/s:\n', vset);
fprintf('%-20s | %-15s | %-40s\n', 'Mode', 'Eigenvalue', 'Dominant States (Scaled Magnitude)');
fprintf('---------------------------------------------------------------------------------------------\n');

names = {'Dutch Roll 1', 'Dutch Roll 2', 'Roll Subsidence', 'Spiral'};
for i = 1:4
    vec = V(:,i);
    % Apply scaling
    scaled_vec = abs(vec .* scaling');
    % Normalize to max=1
    scaled_vec = scaled_vec / max(scaled_vec);
    
    [~, sort_idx] = sort(scaled_vec, 'descend');
    state_names = {'v/V', 'pb/2V', 'rb/2V', 'phi'};
    coupling = sprintf('%s(%.2f), %s(%.2f), %s(%.2f)', ...
        state_names{sort_idx(1)}, scaled_vec(sort_idx(1)), ...
        state_names{sort_idx(2)}, scaled_vec(sort_idx(2)), ...
        state_names{sort_idx(3)}, scaled_vec(sort_idx(3)));
    
    fprintf('%-20s | %-15.4f%+-.4fi | %s\n', names{i}, real(eg(i)), imag(eg(i)), coupling);
end
fprintf('\n');
