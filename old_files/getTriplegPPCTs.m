function [ii,jj,vv1,vv2] = getTriplegPPCTs(CTs,iPs,uoldct1,uoldct2,...
                                    typeCTs,nodeCTs,msh,kk1,kk2,pa,typeG)
% Get triplets of int_CTs(kk*g(uold)*phi*phi)
% Input: - information of NCTs of a specific subdomain
%        - coefficients kk on this subdomain
%        - uold
%        - type of g(u), see file defG.m
% Output: triplet on CTs: vv1 for ij, vv2 for k(i)k(j) (column-arrays)

nCTs = size(CTs,2); % number of all not cut triangles in Omg1
ii = zeros(9*nCTs,1); jj = zeros(9*nCTs,1); 
vv1 = zeros(9*nCTs,1); vv2 = zeros(9*nCTs,1);
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
points = msh.p;

idx = 1;
for t=1:nCTs
    triangle = CTs(:,t); % info of triangle t-th
    uoldT1(1:3) = uoldct1(CTs(1:3,t)); % uold_1's values at vertices
    uoldT2(1:3) = uoldct2(CTs(1:3,t)); % uold_1's values at vertices
    iP1 = iPs(:,1,t); % 1st intersection point
    iP2 = iPs(:,2,t); % 2nd intersection point
    for i=1:3
        for j=1:3
            ii(idx) = CTs(i,t);
            jj(idx) = CTs(j,t);
            if typeCTs(t)==2 % 1 node in Omg2, 2 nodes in Omg1
                rV = points(:,nodeCTs.eachOmg2(1,t)); 
                    % the only node in Omg2
                vWhole1 = getTriplegPPWhole(triangle,uoldT1,j,i,...
                                                    dim,deg,msh,typeG,kk1);
                vPart1 = getTriplegPPPart(triangle,uoldT1,j,i,...
                                    iP1,iP2,rV,dim,deg,typeG,msh,kk1);
                vv1(idx) = vWhole1 - vPart1;
                vv2(idx) = getTriplegPPPart(triangle,uoldT2,j,i,...
                                    iP1,iP2,rV,dim,deg,typeG,msh,kk2);
                idx = idx+1;
            else % typeCT = 0 or 4
                rV = points(:,nodeCTs.eachOmg1(1,t)); 
                    % the only node in Omg1
                vWhole2 = getTriplegPPWhole(triangle,uoldT2,j,i,...
                                                    dim,deg,msh,typeG,kk2);
                vPart2 = getTriplegPPPart(triangle,uoldT2,j,i,...
                                    iP1,iP2,rV,dim,deg,typeG,msh,kk2);
                vv2(idx) = vWhole2 - vPart2;
                vv1(idx) = getTriplegPPPart(triangle,uoldT1,j,i,...
                                    iP1,iP2,rV,dim,deg,typeG,msh,kk1);
                idx = idx+1;
            end
        end
    end
end

end