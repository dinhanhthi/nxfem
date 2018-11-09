function val = getNormL2GfhSTD(msh,err)
% Find semi-norm ||k^{1/2}grad||_L2 for stdFEM
% Copy idea from getNormL2G
% State:
% Input: - triangles tris
%        - error: err 1 x nTs
% Output: scalar value


[i,j,v] = getTripleGGNCTs(msh.t,1,msh);

mA = sparse(i,j,v,msh.nStd,msh.nStd);

val = err*mA*err';
val = sqrt(val);

end