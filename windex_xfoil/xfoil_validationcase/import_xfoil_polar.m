function [alpha, CL, CD, CM] = import_xfoil_polar(filename)
% [alpha, CL, CD, CM] = import_xfoil_polar(filename)
%
% Read an airfoil polar file written by XFoil and return section 
% coefficients over angle of attack (in radian).

    [fid, msg] = fopen(filename, 'r');
    if fid < 0
        error(['Cannot open file: ' filename ' -- ' msg]);
    end
    
    alpha = [];
    CL = [];
    CD = [];
    CM = [];
    
    stage = 1;
    while ~feof(fid)
        
        line = fgetl(fid);
        if stage == 1
            if strfind(line, '------')
                stage = 2;
            end
        else
            [a, n] = sscanf(line, '%f', 7);
            if n == 7
                alpha = [alpha; a(1)*pi/180];
                CL = [CL; a(2)];
                CD = [CD; a(3)];
                CM = [CM; a(5)];
            end
        end
    end
    fclose(fid);
    
end