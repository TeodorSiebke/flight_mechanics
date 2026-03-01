function [Cz, Cdp, Cm, Cza] = lookupfoil(acg, fs, alpha, reynolds)
% [Cz, Cdp, Cm] = lookupfoil(acg, fs, alpha, reynolds)
%
% 

  if nargout < 4
    [Cz, Cdp, Cm] = eval_afgrid(acg.afgrid, alpha, reynolds, fs.delta);
  else
    [Cz, Cdp, Cm, Cza] = eval_afgrid(acg.afgrid, alpha, reynolds, fs.delta);
  end
    
end