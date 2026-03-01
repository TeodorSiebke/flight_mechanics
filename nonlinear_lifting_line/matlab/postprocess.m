function [CL, CD, CM, Fsiv, Fsdp, Ms] = postprocess(acg, fs)
% [CL, CD, CM, Fsiv, Fsdp, Ms] = postprocess(acg, fs)
%
% Compute segment forces, coefficients and moments.
%
% acg : aircraft geometry (struct)
% fs  : flight condition (struct)
% 
% CL, CD, CM : Global lift, drag and pitch moment coefficients
% Fsiv : Segment forces as predicted by the lifting-line model alone
% Fsdp : Segment drag forces obtained from airfoil data
% Ms : Segment moment contributions (Mx, My, Mz) including airfoil data

  % normalization
  S = acg.Sref;
  qoo = 0.5*fs.rho * fs.uref.^2;

  % inviscid force contributions 
  [cxyz, vt] = segment_coefficients(acg, fs, fs.gamma);
  
  lv = acg.pb - acg.pa;
  lseg = sqrt( sum(lv.^2, 2) );
  Svs = acg.chord .* lseg;
  Fsiv = diag(Svs) * qoo * cxyz; 
  
  % viscous drag contributions 
  % alpha = atan2(vt(:,3), vt(:,1));
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
  
  % transform from a/c to wind coordinates 
  % [D, C, L]' = T [Fx Fy Fz]'
  T = body2wind(fs.alpha, fs.beta);
         
  % integrated wind-axes forces and moments 
  Fi = sum(Fs, 1);
  Mi = sum(Ms, 1);
  DCL = T * Fi';
  Mw = T * Mi';
  
  % normalize to get coefficients
  CD = (DCL(1)/qoo + acg.ffsl) / S;  % add fuselage contribution to drag 
  CL = DCL(3) / (qoo * S);
  CM = Mw(2) / (qoo*S*acg.cref);
  
end