function [Cz, Cd, Cm, Cza] = blendplate(spl, alpha)
% [Cz, Cd, Cm, Cza] = blendplate(spl, alpha)
%
% Blend from splined data into coefficients for a separated flat plate
% at high values of angle of attack. This is done to permit the NLL solution 
% procedure to evaluate section coefficients for extremely large AoA, because
% the fsolve() call may sometimes jump there at intermediate points.

    % discard end segments with very steep negative slopes 
    cza_limit = -2.0;
   
    % use spline values when alpha in range
    if is_octave
      spamin = min(spl.breaks);
      spamax = max(spl.breaks);
      
      spd = ppder(spl);
      cx = ppval(spd, spamin); 
      while cx(1) < cza_limit 
        spamin = spamin + 0.1*pi/180;
        cx = ppval(spd, spamin); 
      end 
      cx = ppval(spd, spamax); 
      while cx(1) < cza_limit 
        spamax = spamax - 0.1*pi/180;
        cx = ppval(spd, spamax); 
      end
      %fprintf(1, 'reduced range (%f, %f) to (%f, %f)\n', ...
      %           min(spl.breaks), max(spl.breaks), spamin, spamax); 
   
      i1 = (alpha >= spamin) & (alpha <= spamax);
      sc = ppval(spl, alpha(i1));
      scd = ppval(spd, alpha(i1));
    else 
      spamin = min(spl.knots);
      spamax = max(spl.knots);
      
      spd = fnder(spl);
      cx = fnval(spd, spamin); 
      while cx(1) < cza_limit 
        spamin = spamin + 0.1*pi/180;
        cx = fnval(spd, spamin); 
      end 
      cx = fnval(spd, spamax); 
      while cx(1) < cza_limit 
        spamax = spamax - 0.1*pi/180;
        cx = fnval(spd, spamax); 
      end
      
      i1 = (alpha >= spamin) & (alpha <= spamax);
      sc = fnval(spl, alpha(i1));
      spd = fnder(spl);
      scd = fnval(spd, alpha(i1));
    end
    Cz(i1)  = sc(1,:);
    Cza(i1) = scd(1,:);
    Cd(i1)  = sc(2,:);
    Cm(i1)  = sc(3,:);
    
    % blend into plate coefficients
    if sum(i1) < numel(alpha)
    
        % width of the blending region
        blendrange = pi/9.0;
        if is_octave
          sclo = ppval(spl, spamin);
          schi = ppval(spl, spamax);
          sdlo = ppval(spd, spamin);
          sdhi = ppval(spd, spamax);
        else
          sclo = fnval(spl, spamin);
          schi = fnval(spl, spamax);
          sdlo = fnval(spd, spamin);
          sdhi = fnval(spd, spamax);
        end 
        
        % todo: use cubic hermite spline to stitch segments together
        
        [Czp, Cdp, Cmp, Czap] = separated_plate(alpha);
        i2 = (alpha < spamin);
        [f2, f2da] = perlin_blend(spamin - blendrange, spamin, alpha(i2));
        
        Czl = sclo(1) + sdlo(1) * (alpha(i2) - spamin);
        Cdl = sclo(2) + sdlo(2) * (alpha(i2) - spamin);
        Cml = sclo(3) + sdlo(3) * (alpha(i2) - spamin);
        
        Cz(i2) = (1.0 - f2) .* Czp(i2) + f2 .* Czl;
        Cd(i2) = (1.0 - f2) .* Cdp(i2) + f2 .* Cdl;
        Cm(i2) = (1.0 - f2) .* Cmp(i2) + f2 .* Cml;
        
        Cza(i2) = (0.0 - f2da) .* Czp(i2) + f2da .* Czl ...
                + (1.0 - f2) .* Czap(i2) + f2 .* sdlo(1);
        
        i2 = (alpha > spamax);
        [f2, f2da] = perlin_blend(spamax + blendrange, spamax, alpha(i2));
        
        Czl = schi(1) + sdhi(1) * (alpha(i2) - spamax);
        Cdl = schi(2) + sdhi(2) * (alpha(i2) - spamax);
        Cml = schi(3) + sdhi(3) * (alpha(i2) - spamax);
        
        Cz(i2) = (1.0 - f2) .* Czp(i2) + f2 .* Czl; 
        Cd(i2) = (1.0 - f2) .* Cdp(i2) + f2 .* Cdl;
        Cm(i2) = (1.0 - f2) .* Cmp(i2) + f2 .* Cml;
        
        Cza(i2) = (0.0 - f2da) .* Czp(i2) + f2da .* Czl ...
                + (1.0 - f2) .* Czap(i2) + f2 .* sdhi(1);
        
    end
   
    Cz  = reshape(Cz, size(alpha));
    Cd  = reshape(Cd, size(alpha));
    Cm  = reshape(Cm, size(alpha));
   
end

function [Cz, Cd, Cm, Cza] = separated_plate(alpha)
% [Cz, Cd, Cm] = separated_plate(alpha)
% 
% Evaluate coefficients for a fully separated plate.

    CL = sin(2*alpha);
    Cd = 2*sin(alpha).^2;
    Cz = CL.*cos(alpha) + Cd.*sin(alpha);
    Cm = -Cz * 0.25;
    Cza = 6*sin(alpha).^2 .* cos(alpha) - sin(alpha).*sin(2*alpha) + ...
          2*cos(alpha).*cos(2*alpha);
    
end

function f = saturated_blend(a, b, x)
    f = max(min((x - a)/(b - a), 1.0), 0.0); 
end

function [f,fdx] = perlin_blend(a, b, x)
    c = 1 / (b - a);
    s = max(min((x - a).*c, 1.0), 0.0); 
    f = (3*s.^2 - 2*s.^3); 
    fdx = c*(6*s - 6*s.^2);
    fdx(s <= 0) = 0;
    fdx(s >= 1) = 0;
end