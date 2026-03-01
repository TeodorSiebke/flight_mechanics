function [CL, CD, CC, Cm, Cl, Cn] = coefficients(acg, fs)
% [CL, CD, CC, Cm, Cl, Cn] = coefficients(acg, fs)
% 
% acg : aircraft geometry (struct)
% fs  : flight condition (struct)
% 
% Sum force and moment contributions of all vortex segments
% and normalize to obtain coefficients. All coefficients assume that the 
% geometry definitions in acg are the local design frame (+x backwards, +z 
% upward), then the coefficients are returned compatible with the notation
% in Etkin (+x foward, +z downward).
  
  % body-axes forces and moments 
  [Fsiv, Fsdp, Ms] = forces(acg, fs);
  Fb = sum(Fsiv + Fsdp, 1);
  Mb = sum(Ms, 1);
  
  % to wind axes 
  T = body2wind(fs.alpha, fs.beta);
  Fw = T * Fb';
  Mw = T * Mb';
  
  % normalization
  S = acg.Sref;
  qoo = 0.5*fs.rho * fs.uref.^2;
  
  % lookup fuselage drag 
  if isfield(acg, 'fdragtab');
    logrepm = log( fs.uref / fs.nu );
    ffsl = interp1(acg.lgrtab, acg.fdragtab, logrepm, 'linear', 'extrap');
  else 
    ffsl = acg.ffsl;
  end 
  
  % add fuselage drag contribution 
  CD = Fw(1) / (S*qoo) + ffsl / S;
  CC = Fw(2) / (S*qoo);
  CL = Fw(3) / (S*qoo);
  
  % switch signs for compatibility with Etkin's notation 
  Cl = - Mw(1) / (S*qoo * acg.bref);
  Cm =   Mw(2) / (S*qoo * acg.cref);
  Cn = - Mw(3) / (S*qoo * acg.bref);
  
end 