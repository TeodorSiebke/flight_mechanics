function sps = build_h15()
% sps = build_h15()
%
% Example of how to use assemble_afgrid with computed XFoil polar files.
% This  particular call constructs the interpolation tables for the H15
% airfoil used in the wing lab model.

    dir = '../xfoil/';
    fnames = cell(3, 4);
    reynolds = [3 7 15]' * 1e5;
    delta = [-5 0.0 8 20]' * pi/180;
    fnames{1,1} = [dir 'h15-5re3e5.txt'];
    fnames{2,1} = [dir 'h15-5re7e5.txt'];
    fnames{3,1} = [dir 'h15-5re15e5.txt'];
    
    fnames{1,2} = [dir 'h15re3e5.txt'];
    fnames{2,2} = [dir 'h15re7e5.txt'];
    fnames{3,2} = [dir 'h15re15e5.txt'];
    
    fnames{1,3} = [dir 'h15+8re3e5.txt'];
    fnames{2,3} = [dir 'h15+8re7e5.txt'];
    fnames{3,3} = [dir 'h15+8re15e5.txt'];
    
    fnames{1,4} = [dir 'h15+20re3e5.txt'];
    fnames{2,4} = [dir 'h15+20re7e5.txt'];
    fnames{3,4} = [dir 'h15+20re15e5.txt'];
    
    sps = assemble_afgrid(fnames, reynolds, delta);
    
end