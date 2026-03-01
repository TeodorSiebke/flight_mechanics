% task6_lateral_analysis.m
% Comprehensive Lateral-Directional Stability Analysis for the Windex
% 
% This script analyzes Dutch Roll, Roll Subsidence, and Spiral modes
% across a range of airspeeds and Static Margins.

clear; close all;
addpath('flightsim');
fsm = make_fsim();

% --- PARAMETERS ---
v_sweep = 26:1:80; 
x_np = 0.43;
cref = fsm.cref;
sm_targets = [0.03, 0.09, 0.15];
colors = {'r', 'k', 'b'}; % Aft, Nom, Forward
labels = {'3% SM (Aft)', '9% SM (Nom)', '15% SM (Fwd)'};

all_results = struct();

for s = 1:length(sm_targets)
    current_sm = sm_targets(s);
    xcg = x_np - current_sm * cref;
    fsm.cog(1) = xcg;
    
    fprintf('\n--- Lateral Analysis: SM = %.0f%% (x_cg = %.3f m) ---\n', current_sm*100, xcg);
    speed_data = [];
    
    % Use longitudinal trim as a base (v=0, p=0, r=0, phi=0, psi=0)
    x_lon = [v_sweep(1) 0.4 0.0 0.05 0 1000 0 0.0 0.1]'; 
    prev_eg = [];
    
    for i = 1:length(v_sweep)
        vset = v_sweep(i);
        
        % 1. LONGITUDINAL TRIM
        ivar_lon = [1 2 4 8 9]; ifun_lon = [1 2 3 6 8];
        x_lon(1) = vset;
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
        if norm(ftrim) > 1e-3, continue; end
        
        % 2. 3D STATE EXTRACTION
        % x3d: [u v w p q r phi tet psi xe ye ze de da dr dp]
        x3d = zeros(16,1);
        x3d(1) = x_lon(1); % u
        x3d(3) = x_lon(2); % w
        x3d(5) = x_lon(3); % q
        x3d(8) = x_lon(4); % theta
        x3d(12) = x_lon(6); % alt
        x3d(13) = x_lon(8); % de
        x3d(16) = x_lon(9); % dp
        
        % 3. LATERAL LINEARIZATION
        lat_idx = [2 4 6 7]; % v, p, r, phi
        Alat = zeros(4,4);
        xstep = 1e-5;
        for j = 1:4
            idx = lat_idx(j);
            xh = x3d; xh(idx) = xh(idx) + xstep;
            xmh = x3d; xmh(idx) = xmh(idx) - xstep;
            [xdoth] = fplmod3d(0, xh, fsm);
            [xdotmh] = fplmod3d(0, xmh, fsm);
            Alat(:,j) = (xdoth(lat_idx) - xdotmh(lat_idx)) / (2*xstep);
        end
        
        % 4. EIGENVALUES AND TRACKING
        eg = eig(Alat);
        if isempty(prev_eg)
            % Initial sort: Dutch Roll (complex), Roll (fast real), Spiral (slow real)
            [~, s_idx] = sort(abs(real(eg)), 'descend');
            eg_sorted = eg(s_idx);
        else
            new_idx = zeros(4,1); available = 1:4;
            for k = 1:4
                dist = abs(eg(available) - prev_eg(k));
                [~, best_loc] = min(dist);
                best_idx = available(best_loc);
                new_idx(k) = best_idx; available(best_loc) = [];
            end
            eg_sorted = eg(new_idx);
        end
        prev_eg = eg_sorted;
        
        % Mode analysis for plotting (T1/2 or T2)
        % Note: Dutch Roll is a pair (1,2), Roll is 3, Spiral is 4 if sorted by real part magnitude initially
        t_vals = log(2) ./ real(eg_sorted);
        
        v_idx = length(speed_data) + 1;
        speed_data(v_idx).speed = vset;
        speed_data(v_idx).eig = eg_sorted;
        speed_data(v_idx).t_dr = mean(t_vals(1:2));
        speed_data(v_idx).t_roll = t_vals(3);
        speed_data(v_idx).t_spiral = t_vals(4);
        speed_data(v_idx).freq_dr = mean(abs(imag(eg_sorted(1:2))));
    end
    all_results(s).speed_data = speed_data;
end

% --- PLOTTING ---
figure(3); clf; set(gcf, 'Position', [100 100 1000 500], 'Color', 'w');
figure(4); clf; set(gcf, 'Position', [100 100 700 800], 'Color', 'w');

for s = 1:length(sm_targets)
    sd = all_results(s).speed_data;
    vs = [sd.speed];
    eigs = [sd.eig];
    t_dr = [sd.t_dr];
    t_roll = [sd.t_roll];
    t_spiral = [sd.t_spiral];
    freq_dr = [sd.freq_dr];
    
    % Root Locus
    figure(3);
    subplot(1,2,1); hold on; grid on;
    plot(real(eigs), imag(eigs), [colors{s} '.'], 'HandleVisibility', 'off');
    plot(real(eigs(:,1)), imag(eigs(:,1)), [colors{s} 'o'], 'HandleVisibility', 'off');
    plot(real(eigs(:,end)), imag(eigs(:,end)), [colors{s} 'x'], 'HandleVisibility', 'off');
    plot(nan, nan, colors{s}, 'DisplayName', labels{s});
    
    subplot(1,2,2); hold on; grid on; % Zoom on Spiral/Roll
    plot(real(eigs), imag(eigs), [colors{s} '.'], 'HandleVisibility', 'off');
    plot(nan, nan, colors{s}, 'DisplayName', labels{s});
    axis([-5 0.5 -1 1]);
    
    % Times
    figure(4);
    subplot(3,1,1); hold on; grid on; % Dutch Roll Period or Damping
    plot(vs, t_dr, [colors{s} '-s'], 'DisplayName', labels{s});
    subplot(3,1,2); hold on; grid on; % Roll Subsidence
    plot(vs, t_roll, [colors{s} '-x'], 'DisplayName', labels{s});
    subplot(3,1,3); hold on; grid on; % Spiral
    plot(vs, t_spiral, [colors{s} '-o'], 'DisplayName', labels{s});
end

figure(3);
subplot(1,2,1); title('Lateral Root Locus'); xlabel('Real'); ylabel('Imag');
plot(nan, nan, 'go', 'DisplayName', 'Start: 26 m/s');
plot(nan, nan, 'rx', 'DisplayName', 'End: 80 m/s');
legend('Location', 'best');
subplot(1,2,2); title('Spiral/Roll Region Zoom'); xlabel('Real'); ylabel('Imag');
legend('Location', 'best', 'FontSize', 8);

figure(4);
subplot(3,1,1); title('Dutch Roll: Time to Half ($T_{1/2}$)', 'Interpreter', 'latex'); ylabel('Time (s)');
legend('Location', 'best', 'FontSize', 8);
subplot(3,1,2); title('Roll Subsidence: Time to Half ($T_{1/2}$)', 'Interpreter', 'latex'); ylabel('Time (s)');
legend('Location', 'best', 'FontSize', 8);
subplot(3,1,3); title('Spiral Mode: Growth/Damping ($T_2$ or $T_{1/2}$)', 'Interpreter', 'latex'); xlabel('Airspeed (m/s)'); ylabel('Time (s)');
ylim([-1000 1000]); % Large range for spiral
legend('Location', 'best', 'FontSize', 8);
sgtitle('Lateral-Directional Mode Evolution vs Static Margin');

saveas(3, 'task6_lateral_root_locus.png');
saveas(4, 'task6_lateral_times.png');
fprintf('\nLateral analysis complete.\n');
