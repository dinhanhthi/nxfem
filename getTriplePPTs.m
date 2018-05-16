function [ii,jj,vv] = getTriplePPTs(msh,pa)
% get triple for int_{all tris} phi*phi
% for the level set equation (note 6)
% cf. getGMlsChopp07.m
% NOTE: don't use NXFEM, it's standard FEM
% this function borrows idea from getTriplePPNCTs.m but K is diff
% Input:
% Output: triple ii,jj,vv (column array)

tris = msh.t;
nTs = size(tris,2);
ii = zeros(9*nTs,1); jj = zeros(9*nTs,1); vv = zeros(9*nTs,1);
K=[]; % K is empty and will take the default values in the functions need it

idx=1;
for t=1:nTs
    triangle = tris(:,t);
    for i=1:3
        for j=1:3
            ii(idx) = tris(i,t);
            jj(idx) = tris(j,t);
            vv(idx) = getTriplePPWhole(triangle,j,i,pa,msh,K);
            idx = idx+1;
        end
    end
end

end