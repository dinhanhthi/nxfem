function [ii,jj,vv] = getTripleGGGGTs_gP(vold,msh,del)
% velocity is grad of Phi
% get triples for \int_{all tris} del*(grad v*grad phi)*(grad v*grad phi)
% for the level set equation (note 6)
% cf. getGMlsChopp07.m
% NOTE: don't use NXFEM, it's standard FEM
% Input: - v (already known) = vold
%        - coeff del
% Output: ii,jj,vv

tris = msh.t;
nTs = size(tris,2);
ii = zeros(9*nTs,1); jj = zeros(9*nTs,1); vv = zeros(9*nTs,1);

idx=1;
for t=1:nTs
    for i=1:3
        for j=1:3
            gP1 = getGrad(1,t,msh);
            gP2 = getGrad(2,t,msh);
            gP3 = getGrad(3,t,msh);
            ii(idx) = tris(i,t);
            jj(idx) = tris(j,t);
            gPj = getGrad(j,t,msh);
            gvgPj = vold(tris(1,t))*dot(gP1,gPj)...
                    + vold(tris(2,t))*dot(gP2,gPj)...
                    + vold(tris(3,t))*dot(gP3,gPj); % gvPj is constant
            gPi = getGrad(i,t,msh);
            gvgPi = vold(tris(1,t))*dot(gP1,gPi)...
                    + vold(tris(2,t))*dot(gP2,gPi)...
                    + vold(tris(3,t))*dot(gP3,gPi); % gvPi is constant
            vv(idx) = del*gvgPj*gvgPi;
            idx = idx+1;
        end
    end
end

end
