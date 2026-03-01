function sps = build_finitewing()
% sps = build_finitewing()
%
% Constructs the interpolation tables for the Finite Wing
% airfoils (H15-5).

    % Direct to parent directory of xfoil data
    % Assuming user saves files in a 'polars' dir or same as windex
    dir = '../windex_xfoil/'; 
    
    dir = '../windex_xfoil/'; 
    
    % Dynamically detect available h15_re*e5.txt files
    re_found = [];
    fnames_found = {};
    
    for r = 1:20 % Check from 1e5 to 20e5
        fname = sprintf('%sh15_re%de5.txt', dir, r);
        if exist(fname, 'file') == 2
            re_found = [re_found; r * 1e5];
            fnames_found = [fnames_found; {fname}];
        end
    end

    if ~isempty(re_found)
        reynolds = re_found;
        fnames = fnames_found;
        delta = [ 0.0 ]' * pi/180; % No flaps on this wing currently
        sps = assemble_afgrid(fnames, reynolds, delta);
    else
        % Return dummy if files missing, so main scripts can at least load
        sps = []; 
        warning('No XFOIL polar files found in %s matching h15_re*e5.txt', dir);
    end
    
end
