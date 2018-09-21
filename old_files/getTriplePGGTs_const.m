function [ii,jj,vv] = getTriplePGGTs_const(vold,msh,pa,del)
% velo is a constant
% we use notation vold for velo
% get triples for \int_{all tris} phi*delta*(grad v*grad phi)
% for the level set equation (note 6)
% cf. 
% NOTE: don't use NXFEM, it's standard FEM
% Input: - v (already known) = vold
%        - coeff del
% Output: ii,jj,vv

tris = msh.t;
nTs = size(tris,2);
ii = zeros(9*nTs,1); jj = zeros(9*nTs,1); vv = zeros(9*nTs,1);


idx=1;
for t=1:nTs
    triangle = tris(:,t);
    for i=1:3
        for j=1:3
            ii(idx) = tris(i,t);
            jj(idx) = tris(j,t);
            gPj = getGrad(j,t,msh);
%             gP1 = getGrad(1,t,msh);
%             gP2 = getGrad(2,t,msh);
%             gP3 = getGrad(3,t,msh);
%             gvgPi = vold(tris(1,t))*dot(gP1,gPi)...
%                     + vold(tris(2,t))*dot(gP2,gPi)...
%                     + vold(tris(3,t))*dot(gP3,gPi); % gvPi is constant
            gvgPj = dot(vold,gPj);
            P = []; % force P=1 inside getLWhole
            vv(idx) = del*gvgPj*getLWhole(triangle,i,msh,pa,P);
            idx = idx+1;
        end
    end
end

end