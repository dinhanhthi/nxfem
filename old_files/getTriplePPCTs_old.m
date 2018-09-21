function [ii,jj,vv1,vv2] = getTriplePPCTs(CTs,CT,msh,pa,K)
% Get triplets of int_CTs(K*phi*phi)
% This file is used for general K (the size of K depends on number of
%       Gaussian points (pa.degN)
% Input: - information of NCTs of a specific subdomain
%        - component K1 and K2 for the whole and part case
% Output: triplet on CTs: vv1 for ij, vv2 for k(i)k(j) (column-arrays)

iPs=CT.iPs; typeCTs=CT.type;nodeCTs=CT.nodes;
KW1=K.CTW1; KW2=K.CTW2; KP1=K.CTP1; KP2=K.CTP2;
nCTs = size(CTs,2); % number of all not cut triangles in Omg1
ii = zeros(9*nCTs,1); jj = zeros(9*nCTs,1); 
vv1 = zeros(9*nCTs,1); vv2 = zeros(9*nCTs,1);
points = msh.p;

idx = 1;
for t=1:nCTs
    tri = CTs(:,t); % info of triangle t-th
    iP1 = iPs(:,1,t); % 1st intersection point
    iP2 = iPs(:,2,t); % 2nd intersection point
    KKW1=KW1(t,:); KKW2=KW2(t,:); KKP1=KP1(t,:); KKP2=KP2(t,:);
    for i=1:3
        for j=1:3
            ii(idx) = CTs(i,t);
            jj(idx) = CTs(j,t);
            if typeCTs(t)==2 % 1 node in Omg2, 2 nodes in Omg1
                rV = points(:,nodeCTs.eachOmg2(1,t)); % the only node in Omg2
                vWhole1 = getTriplePPWhole(tri,j,i,pa,msh,KKW1);
                vPart1 = getTriplePPPart(tri,j,i,iP1,iP2,rV,pa,msh,KKP1);
                vv1(idx) = vWhole1 - vPart1;
                vv2(idx) = getTriplePPPart(tri,j,i,iP1,iP2,rV,pa,msh,KKP2);
                idx = idx+1;
            else % typeCT = 0 or 4
                rV = points(:,nodeCTs.eachOmg1(1,t)); % the only node in Omg1
                vWhole2 = getTriplePPWhole(tri,j,i,pa,msh,KKW2);
                vPart2 = getTriplePPPart(tri,j,i,iP1,iP2,rV,pa,msh,KKP2);
                vv2(idx) = vWhole2 - vPart2;
                vv1(idx) = getTriplePPPart(tri,j,i,iP1,iP2,rV,pa,msh,KKP1);
                idx = idx+1;
            end
        end
    end
end

end