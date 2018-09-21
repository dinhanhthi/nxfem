function [ii,jj,vv] = getTripleVeloGPP2(msh,pa,Ts,delT,ri,defV)
% PARALLEL TESTING --> MUCH SLOWER!!!
% Get triples for \int_NCTs del_T*grad_x(xi_i)\int_Ts xi_j*velo_x dx
%   where del_T is optional, grad_x or grad_y of xi_i or xi_j, int of xi_j
%   or xi_i * velo_x or velo_y
% First used for the level set equation (note 6 + hw_levelset_13718.pdf)
% Status: checked with getTriplePGGTs_const.m
% Related: getMEls.m and getMHls.m
% NOTE: don't use NXFEM, it's standard FEM
% Input: - all triangles Ts
%        - del_T: 1 x nTs (SUPG coefficients)
%        - ri = 'xi, yi, xj, yj'
%        - defV(x,y,pa) (function handle)
% Output: ii,jj,vv w.r.t term int_Ts velo\cdot grad phi *graphi 

nTs = size(Ts,2); points = msh.p;
% ii = zeros(9*nTs,1); jj = zeros(9*nTs,1); vv = zeros(9*nTs,1);

gP = getGradPhi(Ts,msh); % grad of phi (2 coordinates x 3 vertices x nTs)

dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points
P = zeros(1,nwt);

for t=1:nTs
    triangle = Ts(:,t);
    v1 = points(:,triangle(1)); % vertex 1
    v2 = points(:,triangle(2)); % vertex 2
    v3 = points(:,triangle(3)); % vertex 3
    for k=1:nwt
        [xk,yk] = getCoorSTD(pt(:,k),v1,v2,v3);
        P(k) = defV(xk,yk,pa);
    end
end

iA = zeros(9*nTs,nTs);
jA = zeros(9*nTs,nTs);
vA = zeros(9*nTs,nTs);
parfor t=1:nTs
    triangle = Ts(:,t);
    del = delT(t); % delta
    it = zeros(9*nTs,1); % column
    jt = zeros(9*nTs,1); % column
    vt = zeros(9*nTs,1); % column
    for i=1:3
        for j=1:3
            idx = (t-1)*9 + (i-1)*3 + j;
            it(idx) = triangle(i);
            jt(idx) = triangle(j);
            switch ri
                case 'xi' % grad_x(xi_i)
                    gxPi = gP(1,i,t);
                    vt(idx) = del*gxPi*getfPhiWhole(msh,pa,triangle,j,P);
                case 'yi'
                    gxPi = gP(2,i,t);
                    vt(idx) = del*gxPi*getfPhiWhole(msh,pa,triangle,j,P);
                case 'xj'
                    gxPi = gP(1,j,t);
                    vt(idx) = del*gxPi*getfPhiWhole(msh,pa,triangle,i,P);
                case 'yj'
                    gxPi = gP(2,j,t);
                    vt(idx) = del*gxPi*getfPhiWhole(msh,pa,triangle,i,P);
            end
        end
    end
    iA(:,t) = it;
    jA(:,t) = jt;
    vA(:,t) = vt;
end

ii = sum(iA,2);
jj = sum(jA,2);
vv = sum(vA,2);

end