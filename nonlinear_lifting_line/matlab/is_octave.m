function r = is_octave()
% r = is_octave()
% 
% Return whether this runs in octave (false in matlab) 

  persistent x;
  if (isempty(x))
    x = exist('OCTAVE_VERSION', 'builtin');
  end
  r = x;
end