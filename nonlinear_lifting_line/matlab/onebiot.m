function [vihat] = onebiot(pa, pb, pc, xh)
% [vihat] = biot(pa, pb, pc, xh)
% 
% Computes velocity induced at pc by a single horseshoe vortex segment 
% with bound segment connecting pa and pb and wake filaments extending along x. 
%
% pa : First vortex point, matrix 1-by-3
% pb : Second vortex point, matrix 1-by-3
% pc : Induction point, matrix n-by-3
% xh : Wake direction vectors, matrix 1-by-3

    ncp = size(pc, 1);
 
    % this works in octave, but not in matlab 
    a = pc - pa;
    b = pc - pb;

    axx = colcross1(a, xh);
    adx = sum( a .* xh, 2 );
    sqa = sum(a.^2, 2);
        
    bxx = colcross1(b, xh);
    bdx = sum(b .* xh, 2);
    sqb = sum(b.^2, 2);
    
    vihat = diag( 1.0 ./ (sqa - sqrt(sqa).*adx) ) * axx  ...
          - diag( 1.0 ./ (sqb - sqrt(sqb).*bdx) ) * bxx;
    
    % add self-induction term only where non-singular
    t1 = sqrt(sqa.*sqb) + sum(a .* b, 2);
    tnz = abs(t1) > eps;
    if sum(tnz) > 0
      axb = colcross(a, b);
      vihat(tnz,:) = vihat(tnz,:) + diag(1.0 ./ t1(tnz)) * axb(tnz,:);    
    end
    
    vihat = vihat / (4*pi);
    
end

function y = colcross1(a, b)
% y = colcross(a, b)
%
% Compute the cross product of a, b with b a single point 

    y = [ (a(:,2).*b(3) - a(:,3).*b(2)) ...
          (a(:,3).*b(1) - a(:,1).*b(3)) ...
          (a(:,1).*b(2) - a(:,2).*b(1)) ];
end

function y = colcross(a, b)
% y = colcross(a, b)
%
% Compute the cross product of a, b assuming that each row of a and b
% is a 3-component vector.
    y = [ (a(:,2).*b(:,3) - a(:,3).*b(:,2)) ...
          (a(:,3).*b(:,1) - a(:,1).*b(:,3)) ...
          (a(:,1).*b(:,2) - a(:,2).*b(:,1)) ];
end

