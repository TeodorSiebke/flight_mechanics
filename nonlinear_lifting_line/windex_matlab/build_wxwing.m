function sps = build_wxwing()
% sps = build_h15()
%
% Example of how to use assemble_afgrid with computed XFoil polar files.
% This  particular call constructs the interpolation tables for the
% Windex main wing without flaps.  (add flaps later)
% Same data is used for delta=0 and delta=1.

    dir = '../windex_xfoil/';
    
    % Reynolds numbers: 3e5, 7e5, 15e5
    reynolds = [3 7 15]' * 1e5;
    
    % Deflections: 0, 10, 20 (degrees)
    delta_deg = [0 10 20]; 
    delta = delta_deg' * pi/180;
    
    n_re = length(reynolds);
    n_delta = length(delta);
    
    fnames = cell(n_re, n_delta);
    
    for i = 1:n_re
        % Reynolds string for filename (e.g., "3e5", "7e5", "15e5")
        if reynolds(i) == 3e5
            re_str = 're3e5';
        elseif reynolds(i) == 7e5
            re_str = 're7e5';
        else
            re_str = 're15e5';
        end
        
        for j = 1:n_delta
            d_val = delta_deg(j);
            
            % Construct filename based on deflection
            if d_val == 0
                % d=0 files are named "wxwing_reXXXX.txt"
                fname = sprintf('wxwing_%s.txt', re_str);
            else
                % Deflected files are named "wxwing_dXX_reXXXX.txt"
                % Handle negative sign in filename (e.g., d-10)
                fname = sprintf('wxwing_d%d_%s.txt', d_val, re_str);
            end
            
            fnames{i, j} = [dir fname];
        end
    end
    
    sps = assemble_afgrid(fnames, reynolds, delta);
end

