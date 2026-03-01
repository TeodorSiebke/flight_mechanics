function [Cz, Cd, Cm, Cza] = eval_afgrid(sps, alpha, reynolds, delta)
% [Cz, Cd, Cm, Cza] = eval_afgrid(sps, alpha, reynolds, delta)
% 
% Evaluate gridded interpolant for a set of angles of attacks,
% Reynolds numbers and flap deflections. This function is typically called 
% by the function which is bound to acg.foil, which in turn needs to 
% account for different airfoils in different spanwise positions.
% eval_afgrid generates data for one particular airfoil type, but can be
% called with many alpha/Reynolds/delta tuples at once. 
%
% sps : struct returned by assemble_afgrid
% alpha : column vector of AoA (radian)
% reynolds : column vector of Reynolds numbers, size must match alpha
% delta : column vector of flap deflections (leave out when not used)
    
    % clamp reynolds number to the ones present in sps, so as to avoid
    % extrapolation in the reynolds number direction  
    logre = 6^2.58 * (log10(reynolds) .^ (-2.58));
    logre = max(min(logre, sps.logre(end)), sps.logre(1));
    logre = reshape(logre, size(alpha));
    
    if nargin < 4
      
      if is_octave
        
        if numel(size(sps.gcz)) > 2
          error('Stored profile data contains flap deflections: Must specify delta vector in eval_afgrid.');
        end
        
        
        Cz = interpn(sps.aip, sps.logre, sps.gcz, alpha, logre);
        if nargout > 1
          Cd = interpn(sps.aip, sps.logre, sps.gcd, alpha, logre);
        end 
        if nargout > 2
          Cm = interpn(sps.aip, sps.logre, sps.gcm, alpha, logre);
        end 
        if nargout > 3
          Cza = interpn(sps.aip, sps.logre, sps.gcza, alpha, logre);
        end 
      else
        Cz = sps.ipcz(alpha, logre);
        Cd = sps.ipcd(alpha, logre);
        Cm = sps.ipcm(alpha, logre);
        Cza = sps.ipcza(alpha, logre);
      end
    else
      if is_octave 
        Cz = interpn(sps.aip, sps.logre, sps.delta, sps.gcz, ... 
                      alpha, logre, delta);
        if nargout > 1
          Cd = interpn(sps.aip, sps.logre, sps.delta, sps.gcd, ... 
                        alpha, logre, delta);
        end 
        if nargout > 2
          Cm = interpn(sps.aip, sps.logre, sps.delta, sps.gcm, ... 
                        alpha, logre, delta);
        end 
        if nargout > 3
          Cza = interpn(sps.aip, sps.logre, sps.delta, sps.gcza, ... 
                        alpha, logre, delta);
        end 
      else
        % Check dimensionality of the interpolant
        if ndims(sps.ipcz.Values) > 2
            Cz = sps.ipcz(alpha, logre, delta);
            Cd = sps.ipcd(alpha, logre, delta);
            Cm = sps.ipcm(alpha, logre, delta);
            Cza = sps.ipcza(alpha, logre, delta);
        else
            Cz = sps.ipcz(alpha, logre);
            Cd = sps.ipcd(alpha, logre);
            Cm = sps.ipcm(alpha, logre);
            Cza = sps.ipcza(alpha, logre);
        end
      end 
    end
    
end