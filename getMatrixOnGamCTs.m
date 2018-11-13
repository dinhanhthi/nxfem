function [ii,jj,vv1,vv2,vv3,vv4]... 
                = getMatrixOnGamCTs(CTs,iPs,uNormalsCT,msh,pa,cp)
% Find all terms on the part of interface that cut triangle
% Input: cut triangles, intersection points
% Output: ii,jj,vv1,vv2,vv3,vv4: triplet vectors to form matrices

nCTs = size(CTs,2); % number of cut triangles
% pa.degP1D: Gaussian quadrature points in 1D (for polynomial function)
dim=1; deg=pa.degP1D; % quadrature's info
kk1 = pa.kk1; kk2 = pa.kk2;
points = msh.p;

% Gradient of function phi_i of the triangles (i is the vertex of triangle)
gradPhiCT = getGradPhi(CTs,msh); % 2 coordinates x 3 vertices x nCTs

ii = zeros(9*nCTs,1); jj = zeros(9*nCTs,1); % column-arrays
vv1 = zeros(9*nCTs,1); vv2 = zeros(9*nCTs,1); % column-arrays
vv3 = zeros(9*nCTs,1); vv4 = zeros(9*nCTs,1); % column-arrays

idx=1;
v = zeros(2,3); % vertices of each triangle
for t=1:nCTs
    pointA = iPs(:,1,t); % 1st intersection point
    pointB = iPs(:,2,t); % 2nd intersection point
    v(:,1) = points(:,CTs(1,t)); % vertx 1
    v(:,2) = points(:,CTs(2,t)); % vertx 2
    v(:,3) = points(:,CTs(3,t)); % vertx 3
    kapp1 = cp.kap1(t); kapp2 = cp.kap2(t);
    lam = cp.lambda(t);
    for i=1:3
        for j=1:3
            ii(idx) = CTs(i,t);
            jj(idx) = CTs(j,t);
            % A_{ij}=a(phi_j,phi_i)
            gradnPhi1 = intGradnPhi(pointA,pointB,gradPhiCT(:,j,t),...
                    uNormalsCT(:,t),i,v,dim,deg);
            gradnPhi2 = intGradnPhi(pointA,pointB,gradPhiCT(:,i,t),...
                    uNormalsCT(:,t),j,v,dim,deg);
            phiphi = intPhiPhi(pointA,pointB,j,i,v,dim,deg);
            vv1(idx) = - kapp1*kk1*gradnPhi1...
                    - kapp1*kk1*gradnPhi2 + lam*phiphi; % A_ij
            vv2(idx) = + kapp2*kk2*gradnPhi1...
                    + kapp2*kk2*gradnPhi2 + lam*phiphi; % A_k(i)k(j)
            vv3(idx) = + kapp1*kk1*gradnPhi1...
                    - kapp2*kk2*gradnPhi2 - lam*phiphi; % A_k(i)j
            vv4(idx) = - kapp2*kk2*gradnPhi1...
                    + kapp1*kk1*gradnPhi2 - lam*phiphi; % A_ik(j)
            idx = idx+1;
        end % end for j
    end % end for i
end % end for t

end