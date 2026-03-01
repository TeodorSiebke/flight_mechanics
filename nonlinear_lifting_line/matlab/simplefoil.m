function [Cz, Cdp, Cm, Cza] = simplefoil(acg, fs, alpha, reynolds)
% [Cz, Cdp, Cm, Cza] = simplefoil(acg, fs, alpha, reynolds)
%
% Extremely primitive airfoil model, which assumes simply inviscid 2D
% lift and moment and returns Reynolds-scaled value of the flat-plate 
% friction drag. There is no pressure drag model, nor does it account for
% flap deflections.
%
% Use this airfoil model only to debug your code! Replace with a proper
% model based on measured data or 2D computations.

    Cz = 2*pi*alpha;
    Cdp = 0.005 * 101.77 * (log10(reynolds) .^ -2.58); 
    Cm = zeros(size(Cz));
    Cza = 2*pi * ones(size(alpha));
    
end