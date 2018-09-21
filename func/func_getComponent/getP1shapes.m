function [rS,drSdr,drSds] = getP1shapes(r,s)
% Define linear shape functions and their partial derivatives at a point
% (r,s) in the "reference triangle"
% Input: coordinate (r,s) in the reference triangle
% Output: the value of basis function and their derivatives at (r,s)

rS = [1-r-s; r; s]; % reference shape functions
drSdr = [-1;1;0]; % partial derivative of rS wrt r
drSds = [-1;0;1]; % partial derivative of rS wrt r
end