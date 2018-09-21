function [ii,jj,vv] = getTripleVeloGPP(Ts,delT,ri,defVelo,dt,pa,msh)
% get triples for \int_NCTs del_T*grad_x(xi_i)\int_Ts xi_j*velo_x dx
%   where del_T is optional, grad_x or grad_y of xi_i or xi_j, int of xi_j
%   or xi_i * velo_x or velo_y
% first used for the level set equation (note 6 + hw_levelset_13718.pdf)
% cf. getMEls.m and getMHls.m
% NOTE: don't use NXFEM, it's standard FEM
% Input: - all triangles Ts
%        - del_T: 1 x nTs (SUPG coefficients)
%        - ri = 'xi, yi, xj, yj'
%        - indicator i or j
%        - velo_x or velo_y (def of velo): depends on x,y,t
% Output: ii,jj,vv w.r.t term int_Omg velo\cdot grad phi *graphi 

nTs = size(Ts,2);
ii = zeros(9*nTs,1); jj = zeros(9*nTs,1); vv = zeros(9*nTs,1);

gP = getGradPhi(Ts,msh); % grad of phi (2 coordinates x 3 vertices x nTs)

    function P = getP(tris,defVelo,dt,pa,msh)
    % get function P used in getLWhole for the case int f(x,y)*Phi
    % Output: matrix P: ntris x nwt
        ntris = size(tris,2);
        points = msh.p;
        dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
        [wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
        nwt = size(wt,2); % number of Gaussian points
        P = zeros(ntris,nwt);
        for tt=1:ntris
            triangle = tris(:,tt);
            v1 = points(:,triangle(1)); % vertex 1
            v2 = points(:,triangle(2)); % vertex 2
            v3 = points(:,triangle(3)); % vertex 3
            for k=1:nwt
                [xk,yk] = getCoorSTD(pt(:,k),v1,v2,v3);
                P(tt,k) = defVelo(xk,yk,dt);
            end
        end
    end

P = getP(Ts,defVelo,dt,pa,msh);

idx=1;
for t=1:nTs
    tri = Ts(:,t);
    del = delT(t); % delta
    PP = P(t,:);
    for i=1:3
        for j=1:3
            ii(idx) = Ts(i,t);
            jj(idx) = Ts(j,t);
            switch ri
                case 'xi' % grad_x(xi_i)
                    gxPi = gP(1,i,t);
                    vv(idx) = del*gxPi*getLWhole(tri,j,msh,pa,PP);
                case 'yi'
                    gxPi = gP(2,i,t);
                    vv(idx) = del*gxPi*getLWhole(tri,j,msh,pa,PP);
                case 'xj'
                    gxPi = gP(1,j,t);
                    vv(idx) = del*gxPi*getLWhole(tri,i,msh,pa,PP);
                case 'yj'
                    gxPi = gP(2,j,t);
                    vv(idx) = del*gxPi*getLWhole(tri,i,msh,pa,PP);
            end
            idx = idx+1;
        end
    end
end

end
