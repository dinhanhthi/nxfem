function [ii,jj,vv] = getTripleGGGGTs(gP,msh,delT)
% Velocity is grad of Phi
% Get triples for \int_{all tris} del*(grad v*grad phi)*(grad v*grad phi)
% For the level set equation (note 6)
% cf. main_chopp*, getMElsgP, getMHlsgP
% NOTE: don't use NXFEM, it's standard FEM
% Input: - gP: gP.x (1 x nTs), gP.y (1 x nTs): grad of v
%        - coeff del
% Output: ii,jj,vv
% ----------------------------------------------------------
% Update 26/10/18: The OLD file used to compute grad v on vertices of each
%   triangle. We now use pdegrad to find grad v on the center of each
%   triangle (this file), i.e. grad v is constant on whole triangle and discont at edge.
% ----------------------------------------------------------

tris = msh.t;
nTs = size(tris,2);
ii = zeros(9*nTs,1); jj = zeros(9*nTs,1); vv = zeros(9*nTs,1);

gradPhi = getGradPhi(tris,msh); % 2 coor x 3 vertices x nTris

idx=1;
for t=1:nTs
    del = delT(t); % delta
    gv = [gP.x(t), gP.y(t)];
    for i=1:3
        for j=1:3
            ii(idx) = tris(i,t);
            jj(idx) = tris(j,t);
            gvgPj = dot(gv,gradPhi(:,j,t)); % gvPj is constant
            gvgPi = dot(gv,gradPhi(:,i,t)); % gvPi is constant
            vv(idx) = del*gvgPj*gvgPi;
            idx = idx+1;
        end
    end
end

end
