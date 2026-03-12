function sps = build_wxtail()
% sps = build_h15()
%
% Example of how to use assemble_afgrid with computed XFoil polar files.
% This  particular call constructs the interpolation tables for the Windex
% tail plane airfoils.  Add more control surface deflections to improve.

    dir = '../windex_xfoil/';
    fnames = cell(3, 4);
    reynolds = [3 7 15]' * 1e5;
    delta = [-20.0 -10.0 0.0 10.0]' * pi/180;
    
    fnames{1,1} = [dir 'wxtail_d-20_re3e5.txt'];
    fnames{2,1} = [dir 'wxtail_d-20_re7e5.txt'];
    fnames{3,1} = [dir 'wxtail_d-20_re15e5.txt'];

    fnames{1,2} = [dir 'wxtail-10re3e5.txt'];
    fnames{2,2} = [dir 'wxtail-10re7e5.txt'];
    fnames{3,2} = [dir 'wxtail-10re15e5.txt'];
    
    fnames{1,3} = [dir 'wxtail+0re3e5.txt'];
    fnames{2,3} = [dir 'wxtail+0re7e5.txt'];
    fnames{3,3} = [dir 'wxtail+0re15e5.txt'];
    
    fnames{1,4} = [dir 'wxtail+10re3e5.txt'];
    fnames{2,4} = [dir 'wxtail+10re7e5.txt'];
    fnames{3,4} = [dir 'wxtail+10re15e5.txt'];
    
%    fnames{1,4} = [dir 'h15+20re3e5.txt'];
%    fnames{2,4} = [dir 'h15+20re7e5.txt'];
%    fnames{3,4} = [dir 'h15+20re15e5.txt'];
    
    sps = assemble_afgrid(fnames, reynolds, delta);
    
end
