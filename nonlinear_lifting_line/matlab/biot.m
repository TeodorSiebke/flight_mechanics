function [vihat] = biot(pa, pb, pc, xh)
% [vihat] = biot(pa, pb, pc, xh)
% 
% Computes velocity induced at pc by a horseshoe vortex with bound segment
% connecting pa and pb and wake filaments extending along x. 
%
% pa : First vortex point, matrix n-by-3
% pb : Second vortex point, matrix n-by-3
% pc : Induction point, matrix n-by-3
% xh : Wake direction vectors, matrix n-by-3

    nvx = size(pa, 1);
    ncp = size(pc, 1);
 
    a = pc - pa;
    b = pc - pb;

    axx = colcross(a, xh);
    adx = sum( a .* xh, 2 );
    sqa = sum(a.^2, 2);
        
    bxx = colcross(b, xh);
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

function y = colcross(a, b)
% y = colcross(a, b)
%
% Compute the cross product of a, b assuming that each row of a and b
% is a 3-component vector.
    y = [ (a(:,2).*b(:,3) - a(:,3).*b(:,2)) ...
          (a(:,3).*b(:,1) - a(:,1).*b(:,3)) ...
          (a(:,1).*b(:,2) - a(:,2).*b(:,1)) ];
end

%function [xx,adx,asq] = cc3(pa, pc, xh)
%% [xx,adx,asq] = cc3(pa, pc, xh)
%%
%% Equivalent to
%% a  = pc - repmat(pa, n, 1);
%% xn = repmat(xh, n, 1); 
%% xx = colcross(a, xh);
%% adx = sum(a.*xh, 2);
%% asq = sum(a.^2, 2);
%% 
%    xx = [ ((pc(:,2) - pa(2))*xh(3) - (pc(:,3) - pa(3))*xh(2)) ...
%          ((pc(:,3) - pa(3))*xh(1) - (pc(:,1) - pa(1))*xh(3)) ...
%          ((pc(:,1) - pa(1))*xh(2) - (pc(:,2) - pa(2))*xh(1)) ];
%    adx = (pc(:,1) - pa(1))*xh(1) + (pc(:,2) - pa(2))*xh(2) + (pc(:,3) - pa(3))*xh(3);
%    asq = (pc(:,1) - pa(1)).^2 + (pc(:,2) - pa(2)).^2 + (pc(:,3) - pa(3)).^2;
%end