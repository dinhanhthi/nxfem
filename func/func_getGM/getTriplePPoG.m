function [ii,jj,vv1,vv2,vv3,vv4]= getTriplePPoG(CTs,iPs,msh,pa,L)
% Find triplet for + int_G L*phi*phi (on interface)
% L: 4 x nCTs
% Input: - cut triangles, intersection points
%        - diff coef wrt each triangle 
%        - pa: contains all basis parameters
%        - L: 4 x nCTs
% Output: ii,jj,vv1,vv2,vv3,vv4: triplet vectors to form matrices
%           (column array)

nCTs = size(CTs,2); % number of cut triangles
points = msh.p;

ii = zeros(9*nCTs,1); jj = zeros(9*nCTs,1); % column-arrays
vv1 = zeros(9*nCTs,1); vv2 = zeros(9*nCTs,1); % column-arrays
vv3 = zeros(9*nCTs,1); vv4 = zeros(9*nCTs,1); % column-arrays

if isempty(L)
    L = ones(4,nCTs);
end

idx=1;
v = zeros(2,3); % vertices of each triangle
for t=1:nCTs
    pointA = iPs(:,1,t); % 1st intersection point
    pointB = iPs(:,2,t); % 2nd intersection point
    v(:,1) = points(:,CTs(1,t)); % vertx 1
    v(:,2) = points(:,CTs(2,t)); % vertx 2
    v(:,3) = points(:,CTs(3,t)); % vertx 3
    for i=1:3
        for j=1:3
            ii(idx) = CTs(i,t);
            jj(idx) = CTs(j,t);
            % A_{ij}=a(phi_j,phi_i)
            phiphi = intPhiPhi(pointA,pointB,j,i,v,pa);
            vv1(idx) = L(1,t)*phiphi; % A_ij
            vv2(idx) = L(2,t)*phiphi; % A_k(i)k(j)
            vv3(idx) = - L(3,t)*phiphi; % A_k(i)j
            vv4(idx) = - L(4,t)*phiphi; % A_ik(j)
            idx = idx+1;
        end % end for j
    end % end for i
end % end for t

end