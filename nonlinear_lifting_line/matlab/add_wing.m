function [acg,idx] = add_wing(acg, p1, p2, c1, c2, twist1, twist2, n)
% [acg,idx] = add_wing(acg, p1, p2, c1, c2, n)
% 
% Add an additional wing panel to acg. Chord and twist will be interpolated
% linearly in th espanwise direction. To create a more complex planform,
% call this function repeatedly or write your own version which generates
% vortex segment endpoints in a similar manner. 
%
% p1 : Point at leading edge left tip (1x3)
% p2 : Point at leading edge right tip (1x3)
% c1 : Chord at left tip
% c2 : Chord at right tip 
% twist1 : twist angle at left tip (radian)
% twist2 : twist anegle at right tip (radian)
% n  : number of vortex segments to generate 

  np = n+1;
  eta = linspace(0.0, 1.0, np)';
  chord = (1-eta)*c1 + eta*c2;
  twist = (1-eta)*twist1 + eta*twist2;
  pts = zeros(np, 3);
  for ki = 1:3
    pts(:,ki) = (1-eta)*p1(ki) + eta*p2(ki);
  end 
  [acg,idx] = add_segments(acg, pts, chord, twist);
  
end 