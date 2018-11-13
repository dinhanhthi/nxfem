function [rS,drSdr,drSds] = getP2shapes(r,s)
% Define quadratic shape functions and their partial derivatives at a point
% (r,s) in the "reference triangle"
% Input: coordinate (r,s) in the reference triangle
% Output: the value of basis function and their derivatives at (r,s)

% reference shape functions in P2
rS = [1-3*r-3*s+2*r^2+4*r*s+2*s^2;
      2*r^2-r;
      2*s^2-s;
      4*r*s;
      4*s-4*r*s-4*s^2;
      4*r-4*r^2-4*r*s];
  
% partial derivative of rS wrt r
drSdr = [-3+4*r+4*s; 4*r-1; 0; 4*s; -4*s; 4-8*r-4*s];

% partial derivative of rS wrt r
drSds = [-3+4*r+4*s; 0; 4*s-1; 4*r; 4-4*r-8*s; -4*r]; 
end