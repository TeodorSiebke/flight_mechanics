function [cxyz, vt, Jcz, Jvz] = segment_coefficients(acg, fs, Gamma)
% [cxyz, vt] = segment_coefficients(acg, fs, Gamma)
% 
% Compute force coefficients and effective local angles of attack at segment 
% centers. 
% 
% fs.uref : scalar reference velocity for normalization 
% acg.chord : wing chord at segment centers (ns,1)
% acg.pa, pb : left/right end points of vortex segments (ns,3)
% acg.xh : wake vortex filament directions (ns,3)
% fs.vb : body velocity (ac motion) at segment centers (ns,3)
% Gamma : vertex strengths (ns,1)
%
% cxyz : force coefficients Cx, Cy, Cz in body coordinates (ns,3)
% vt : effective local velocities at segment midpoints (ns,3)

    % number of vortex segments
    nvx = size(acg.pa, 1);
    
    % determine induced velocities at segment midpoints
    pmp = 0.5*(acg.pa + acg.pb);
    vi = zeros(nvx, 3);
    
    % Jacobian of the induced velocity in z-direction 
    % For all Jacobian terms: rows of J correspond to collocation point, 
    % columns of J to the Gamma-values of the vortices, so 
    % dvx(ki, kj) is x-velocity induced at collocation point ki by a vortex
    % segment kj with Gamma(kj) = 1   
    dvx = zeros(nvx, nvx);
    dvy = zeros(nvx, nvx);
    dvz = zeros(nvx, nvx);
    
    % sum contributions to induced velocity vi from all vortex segments ki
    lv = acg.pb - acg.pa;
    Sqoo = 0.5*fs.uref^2 .* acg.chord .* sqrt( sum(lv.^2, 2) );
    
    if is_octave
      persistent have_cbiot;
      if (isempty(have_cbiot))
        have_cbiot = exist('cbiot.oct', 'file');
      end
      if have_cbiot ~= 0
        v3 = cbiot(acg.pa, acg.pb, pmp, acg.xh);
        dvx = v3(:,1:3:end)';
        dvy = v3(:,2:3:end)';
        dvz = v3(:,3:3:end)';
        vi = reshape(Gamma' * v3, [3,nvx])';
      else 
        for ki = 1:nvx
          v = onebiot(acg.pa(ki,:), acg.pb(ki,:), pmp, acg.xh(ki,:));
          vi = vi + Gamma(ki) * v;
          dvx(:,ki) = v(:,1); 
          dvy(:,ki) = v(:,2); 
          dvz(:,ki) = v(:,3); 
        end
      end
    else
      for ki = 1:nvx
        v = biot( repmat(acg.pa(ki,:), nvx, 1), ...
                  repmat(acg.pb(ki,:), nvx, 1), ...
                  pmp, repmat(acg.xh(ki,:), nvx, 1) );
        vi = vi + Gamma(ki) * v;
        dvx(:,ki) = v(:,1); 
        dvy(:,ki) = v(:,2); 
        dvz(:,ki) = v(:,3); 
      end
    end
    
    % effective total velocity 
    vt = vi - fs.vb;
    
    % Jacobian of the velocity component that is "upward" with respect to Gamma
    % For vortex elements parallel to y and with the convention xh = (1,0,0),
    % this will extract Jvz = dvz, for collocation points on the VTP, J = dvy.
    Jvz = diag(acg.zup(:,1)) * dvx + diag(acg.zup(:,2)) * dvy + ...
          diag(acg.zup(:,3)) * dvz;
    
    % segment force coefficient (nvx x 3)
    cxyz = diag(1.0 ./ Sqoo) * cross(vt, diag(Gamma)*lv);
    
    vx = vt(:,1);
    vy = vt(:,2);
    vz = vt(:,3);
    lix = lv(:,1);
    liy = lv(:,2);
    liz = lv(:,3);    
    
    % Jacobian of z-force coefficients
    % 
    % [vx vy vz] x [lix liy liz] * Gi
    % cx = (vy*liz - vz*liy) * Gi
    % cy = (vz*lix - vx*liz) * Gi
    % cz = (vx*liy - vy*lix) * Gi
    %
    % dcx = (vy*liz - vz*liy) + (dvy*liz - dvz*liy) * Gi
    % dcy = (vz*lix - vx*liz) + (dvz*lix - dvx*liz) * Gi
    % dcz = (vx*liy - vy*lix) + (dvx*liy - dvy*lix) * Gi
         
    Jfx = diag(vy.*liz - vz.*liy) + (dvy*diag(liz) - dvz*diag(liy))*diag(Gamma);
    Jfy = diag(vz.*lix - vx.*liz) + (dvz*diag(lix) - dvx*diag(liz))*diag(Gamma);
    Jfz = diag(vx.*liy - vy.*lix) + (dvx*diag(liy) - dvy*diag(lix))*diag(Gamma);
          
    % extract the projection onto the local up-direction and normalize
    Jcz = diag(1.0 ./ Sqoo) * (diag(acg.zup(:,1)) * Jfx + ...
                               diag(acg.zup(:,2)) * Jfy + ...
                               diag(acg.zup(:,3)) * Jfz);
                               
    %    Jcz = diag(1.0 ./ Sqoo) * ...
    %          diag(vt(:,1).*lv(:,2) - vt(:,2).*lv(:,1)) + ...
    %          (dvx*diag(lv(:,2)) - dvy*diag(lv(:,1))) * diag(Gamma);
    
end  