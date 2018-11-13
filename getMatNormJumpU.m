function mNormJumpU = getMatNormJumpU(CTs,iPs,msh,pa,coef)
% This function is to compute the matrix of finding value of norm of [uh] or
% [u-uh] on Gamma. It's first used in finding the value of zeta_S and norm_h,Gam
% given on page 33, Barrau's thesis
% Note that the algorithm studies from function getMatrixOnGamCTs
% Input: - cut triangles CTs
%        - intersection points iPs
%        - coefficients coef, it depends on the triangles (eg. getLambda)
% Output: - a matrix


nCTs = size(CTs,2); % number of cut triangles
% pa.degP1D: Gaussian quadrature points in 1D (for polynomial function)
dim=1; deg=pa.degP1D; % quadrature's info
points=msh.p; newNodes=msh.newNodes;

% setting up
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
    for i=1:3
        for j=1:3
            ii(idx) = CTs(i,t);
            jj(idx) = CTs(j,t);
            % A_{ij}=a(phi_j,phi_i)
            phiphi = intPhiPhi(pointA,pointB,j,i,v,dim,deg);
            vv1(idx) = coef(t)*phiphi; % ij
            vv2(idx) = coef(t)*phiphi; % k(i)k(j)
            vv3(idx) = -coef(t)*phiphi; % k(i)j
            vv4(idx) = -coef(t)*phiphi; % ik(j)
            idx = idx+1;
        end % end for j
    end % end for i
end % end for t

% ij
in=ii; jn=jj; vn=vv1;
% k(i)k(j)
itmp = newNodes(ii); jtmp = newNodes(jj);
in=[in;itmp]; jn=[jn;jtmp]; vn=[vn;vv2]; % column-array
% k(i)j
it3 = newNodes(ii); % column-array
in = [in;it3]; jn = [jn;jj]; vn = [vn;vv3]; % column-array
% ik(j)
jt4 = newNodes(jj); % column-array
in = [in;ii]; jn = [jn;jt4]; vn = [vn;vv4]; % column-array

mNormJumpU = sparse(in,jn,vn);
end