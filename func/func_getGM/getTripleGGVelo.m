function [ii,jj,vv] = getTripleGGVelo(msh,pa,Ts,delT,gi,gj,defV)
% Get triples for \int_NCTs del_T*grad_x(xi_j)*grad_x(xi_i)\int_Ts defV dx
%   where del_T is optional, grad_x or grad_y of xi_i or xi_j
% Status: checked with del=0, 
% First used for the level set equation (note 6 + hw_levelset_13718.pdf)
% cf. getMEls.m and getMHls.m
% NOTE: don't use NXFEM, it's standard FEM
% Input: - all triangles Ts
%        - del_T: 1 x nTs (SUPG coefficients)
%        - ri = 'xi, yi, xj, yj'
%        - defV(x,y,pa) (function handle)
% Output: ii,jj,vv 

nTs = size(Ts,2); points = msh.p;
ii = zeros(9*nTs,1); jj = zeros(9*nTs,1); vv = zeros(9*nTs,1);
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points
gP = getGradPhi(Ts,msh); % grad of phi (2 coordinates x 3 vertices x nTs)


idx=1;
for t=1:nTs
    triangle = Ts(:,t);
    del = delT(t); % delta
    v1 = points(:,triangle(1)); % vertex 1
    v2 = points(:,triangle(2)); % vertex 2
    v3 = points(:,triangle(3)); % vertex 3
    areaT = getAreaTri(v1,v2,v3);
    val = 0;
    for k=1:nwt
        [xk,yk] = getCoorSTD(pt(:,k),v1,v2,v3);
        val = val + defV(xk,yk,pa)*areaT*wt(k);
    end
    for i=1:3
        for j=1:3
            ii(idx) = Ts(i,t);
            jj(idx) = Ts(j,t);
            switch gi
                case 'xi' % grad_x(xi_i)
                    gPi = gP(1,i,t);
                case 'yi'
                    gPi = gP(2,i,t);
            end
            switch gj
                case 'xj' % grad_x(xi_j)
                    gPj = gP(1,j,t);
                case 'yj'
                    gPj = gP(2,j,t);
            end
            vv(idx) = del*gPi*gPj*val;
            idx = idx+1;
        end
    end
end

end