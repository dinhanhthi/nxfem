function [ii,jj,vv1,vv2,vv3,vv4] = getMatEstEtaSonGam(CTs,iPs,msh,pa)
% Get triplet on interface segments for finding matrix etaS
% Input: given in the input parameters
% Output: triplet to build the final matrix

nCTs = size(CTs,2); % number of cut triangles
points = msh.p; % points of the mesh
kk1 = pa.kk1; kk2 = pa.kk2;
kmin = min(kk1,kk2);

ii = zeros(9*nCTs,1); jj = zeros(9*nCTs,1); % column-arrays
vv1 = zeros(9*nCTs,1); vv2 = zeros(9*nCTs,1); % column-arrays
vv3 = zeros(9*nCTs,1); vv4 = zeros(9*nCTs,1); % column-arrays

idx=1;
v = zeros(2,3); % vertices of each triangle
for t=1:nCTs
    mt = CTs(5,t); % idx of triangle t in msh.t
    pointA = iPs(:,1,t); % 1st intersection point
    pointB = iPs(:,2,t); % 2nd intersection point
    % length of edge
    lenEdge = sqrt((pointA(1)-pointB(1))^2 +(pointA(2)-pointB(2))^2);
    v(:,1) = points(:,CTs(1,t)); % vertx 1
    v(:,2) = points(:,CTs(2,t)); % vertx 2
    v(:,3) = points(:,CTs(3,t)); % vertx 3
    % 2 intersection points in the reference triangle
    [xHa,yHa] = getCoorRef(pointA,v(:,1),v(:,2),v(:,3));
    [xHb,yHb] = getCoorRef(pointB,v(:,1),v(:,2),v(:,3)); 
    % length of edge in ref triangle
    lenEdgeRef = sqrt((xHa-xHb)^2 +(yHa-yHb)^2);
    for i=1:3
        for j=1:3
            ii(idx) = CTs(i,t);
            jj(idx) = CTs(j,t);
            % (grad_ni * grad_nj)
            gNgN = getGradnPhi(i,mt,pointA,pointB,msh)*...
                            getGradnPhi(j,mt,pointA,pointB,msh);
            % A_{ij}=a(phi_j,phi_i)
            vv1(idx) = lenEdge*kk1^2/(kmin*lenEdgeRef)*gNgN; % ij
            vv2(idx) = lenEdge*kk2^2/(kmin*lenEdgeRef)*gNgN; % k(i)k(j)
            vv3(idx) = -lenEdge*kk1*kk2/(kmin*lenEdgeRef)*gNgN; % k(i)j
            vv4(idx) = -lenEdge*kk2*kk1/(kmin*lenEdgeRef)*gNgN; % ik(j)
            idx = idx+1;
        end % end for j
    end % end for i
end % end for t
end