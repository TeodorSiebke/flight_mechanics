function [acg, idx] = add_segments(acg, pts, chord, twist)
% [acg, idx] = add_segments(acg, pts, chord, twist)
% 
% acg: struct containing aircraft description (can be empty)
% pts: n-by-3 matrix of points along the leading edge 
% chord: n-vector of chord lengths
% twist: n-vector of twist angles (leave out for zero twist)
%
% idx: indices of the newly created vortex segments in acg (for plotting)

  np = size(pts,1);
  ns = np-1;
  if size(pts,2) ~= 3
    error('pts must be n-by-3 matrix containing points.');
  end 
  if numel(chord) ~= np 
    error('Array of chord values does not match number of points.');
  end 
  if nargin < 4
    twist = zeros(size(chord));
  end 
  
  % generate segments from points 
  ra = 1:ns;
  rb = 2:np;
  xh = repmat([1 0 0], [ns 1]);
  cseg = 0.5*(chord(ra) + chord(rb));
  tseg = 0.5*(twist(ra) + twist(rb));
  
  pa = pts(ra,:) + diag(0.25*chord(ra)) * xh;
  pb = pts(rb,:) + diag(0.25*chord(rb)) * xh;
  
  if isfield(acg, 'chord')
    ibegin = numel(acg.chord) + 1;
    acg.chord = [acg.chord; cseg];
    acg.incidence = [acg.incidence; tseg];
    acg.pa = [acg.pa; pa];
    acg.pb = [acg.pb; pb];
    acg.xh = [acg.xh; xh];
  else
    ibegin = 1; 
    acg.chord = cseg;
    acg.incidence = tseg;
    acg.pa = pa;
    acg.pb = pb;
    acg.xh = xh;
  end 
  idx = ibegin:(ibegin+ns-1);
  
  % compute local 'up' direction (+z) for all vortex segments 
  acg.zup = cross(acg.pa - acg.pb, acg.xh);
  len = sqrt( sum( acg.zup.^2, 2 ) ); 
  acg.zup = diag(1.0 ./ len) * acg.zup;
  
end 