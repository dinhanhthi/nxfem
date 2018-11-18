function [ii,jj,vv1,vv2,vv3,vv4] = getTripleGPoG(CTs,iPs,uNCTs,msh,pa,cp)
% Find triplet for 2 terms int_G k*grad*phi (on interface)
% - int_Gam gradnphi - int_Gam granphi
% Input: - cut triangles, intersection points
%        - diff coef wrt each triangle 
%        - pa: contains all basis parameters
%        - cp: contains parameters wrt the model (kapi,kki)
% Output: ii,jj,vv1,vv2,vv3,vv4: triplet vectors corresponding to 4 cases
%           of phi_i and phi_j (column array)

nCTs = size(CTs,2); % number of cut triangles
kk1 = cp.kk1; kk2 = cp.kk2; 
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
    for i=1:3
        for j=1:3
            ii(idx) = CTs(i,t);
            jj(idx) = CTs(j,t);
            % A_{ij}=a(phi_j,phi_i)
            gradnPhi1 = intGradnPhi(pointA,pointB,gradPhiCT(:,j,t),uNCTs(:,t),i,v,pa); % ji
            gradnPhi2 = intGradnPhi(pointA,pointB,gradPhiCT(:,i,t),uNCTs(:,t),j,v,pa); % ij
            vv1(idx) = -kapp1*kk1*gradnPhi1-kapp1*kk1*gradnPhi2; % A_ij
            vv2(idx) = +kapp2*kk2*gradnPhi1+kapp2*kk2*gradnPhi2;% A_k(i)k(j)
            vv3(idx) = +kapp1*kk1*gradnPhi1-kapp2*kk2*gradnPhi2; % A_k(i)j
            vv4(idx) = -kapp2*kk2*gradnPhi1+kapp1*kk1*gradnPhi2; % A_ik(j)
            idx = idx+1;
        end % end for j
    end % end for i
end % end for t

end