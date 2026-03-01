function T = body2wind(alpha, beta)
  Ta = [ cos(alpha) 0            sin(alpha) 
         0             1            0  
        -sin(alpha) 0            cos(alpha) ];
  Tb = [ cos(beta) -sin(beta) 0 
         sin(beta)  cos(beta) 0 
         0             0            1 ];
  T = Ta * Tb;
end 