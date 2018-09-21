function [mNCTs1,mNCTs2] = getMatrixNCTs(omg1NCTs,omg2NCTs,msh,pa)
% Get matrix on not cut triangles
% Input: triangles on each subdomain
% Output: 2 matrices wrt each subdomain.

[i1,j1,v1] = getMatrixGradGradNCTs(omg1NCTs,pa.kk1,msh); % NCTs1
[i2,j2,v2] = getMatrixGradGradNCTs(omg2NCTs,pa.kk2,msh); % NCTs2

mNCTs1 = sparse(i1,j1,v1); % in NCTs1
mNCTs2 = sparse(i2,j2,v2); % in NCTs2
end