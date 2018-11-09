function val = getNormL2std(msh,pa,err)
% Get the ||err_h||_L2 of a stdFEM solution
% Note: it's diff from getNormL2nxfem which is for nxfem solution
% State: - change fater rewrite to new getTriplePPNCTs & getTriplePPCTs
%        - checked with the old method getMatrixL2.m
% Input: - err in stdFEM: 1 x nTs
%        - mesh info: include all triangles msh.t
% Output: value of L2 error


K = []; % K takes default value, 1
[i,j,v] = getTriplePPNCTs(msh,pa,msh.t,K);
mA = sparse(i,j,v);


val = err*mA*err';
val = sqrt(val);


end