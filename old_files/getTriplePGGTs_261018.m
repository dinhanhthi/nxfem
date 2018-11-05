function [ii,jj,vv] = getTriplePGGTs(vold,msh,pa,delT)
% Get triples for \int_{all tris} phi*delta*(grad v*grad phi)
% Velocity is grad of Phi
% for the level set equation (note 6)
% cf. main_chopp2007, getMEls_gP, getMHls_gP
% REMARK: don't use NXFEM, it's standard FEM
% Input: - v (already known) = vold in STD
%        - del_T: 1 x nTs (SUPG coefficients)
% Output: ii,jj,vv
% ----------------------------------------------------------
% Update 26/10/18: This file used to compute grad v on vertices of each
%   triangle. We now use pdegrad to find grad v on the center of each
%   triangle, i.e. grad v is constant on whole triangle and discont at edge.
% ----------------------------------------------------------

tris = msh.t;
nTs = size(tris,2);
ii = zeros(9*nTs,1); jj = zeros(9*nTs,1); vv = zeros(9*nTs,1);


idx=1;
for t=1:nTs
    triangle = tris(:,t);
    del = delT(t); % delta
    for i=1:3
        for j=1:3
            ii(idx) = tris(i,t);
            jj(idx) = tris(j,t);
            gPi = getGrad(i,t,msh);
            gP1 = getGrad(1,t,msh);
            gP2 = getGrad(2,t,msh);
            gP3 = getGrad(3,t,msh);
            gvgPi = vold(tris(1,t))*dot(gP1,gPi)...
                    + vold(tris(2,t))*dot(gP2,gPi)...
                    + vold(tris(3,t))*dot(gP3,gPi); % gvPi is constant
            P = []; % force P=1 inside getfPhiWhole
            vv(idx) = del*gvgPi*getfPhiWhole(msh,pa,triangle,j,P);
            idx = idx+1;
        end
    end
end

end