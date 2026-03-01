function [Fsiv, Fsdp, Ms] = forces(acg, fs)
% [Fsiv, Fsdp, Ms] = postprocess(acg, fs)
%
% Compute segment forces and moments.
%
% acg : aircraft geometry (struct)
% fs  : flight condition (struct)
% 
% Fsiv : Segment forces as predicted by the lifting-line model alone
% Fsdp : Segment drag forces obtained from airfoil data
% Ms : Segment moment contributions (Mx, My, Mz) including airfoil data

  % inviscid force contributions 
  [cxyz, vt] = segment_coefficients(acg, fs, fs.gamma);
  
  % segment dimensions
  lv = acg.pb - acg.pa;
  lseg = sqrt( sum(lv.^2, 2) );
  Svs = acg.chord .* lseg;
  
  qoo = 0.5*fs.rho * fs.uref.^2;
  Fsiv = diag(Svs) * qoo * cxyz; 
  
  % viscous drag contributions 
  alpha = vortex_alpha(acg, vt);
  
  uc = sqrt(sum( vt.^2, 2 ));
  Re = uc .* acg.chord ./ fs.nu;
  [~, Cdp, Cm] = feval(acg.foil, acg, fs, alpha, Re);
  
  % directed drag forces (aligned with local flow direction!)
  vsn = sqrt(sum(vt.^2, 2));
  Fsdp = diag(Cdp .* Svs .* 0.5 * fs.rho .* vsn) * vt;
  Fs = Fsiv + Fsdp;
  
  % moments about reference point, in vortex coordinates
  nvx = size(acg.pa, 1);  
  r = 0.5*(acg.pa + acg.pb) - repmat(acg.pref, [nvx, 1]);
  Ms = cross(r, Fs);
  
  % add contribution of section pitch moments, which are directed along
  % their respective spanwise axis (pb - pa); these contributions would 
  % be zero for symmetric airfoils in inviscid flow.
  Mpy = Cm .* Svs .* 0.5*fs.rho .* sum(vt.^2, 2);
  Ms = Ms + diag(Mpy ./ lseg) * lv; 
  
end