function [ii,jj,vv1,vv2] = getTriplePPCTs(msh,pa,CTs,CT,K)
% Get triplets of int_CTs(K*phi*phi)
% This file is used for general K (the size of K depends on number of
%       Gaussian points (pa.degN)
% Note: rewriten after 1st rewrite with inputParser (worst performance)
% State: - checked with getGMgPP
% Input: - cut triangles CTs
%        - CTs' info CT
%        - K.CTW1,K.CTW2 for the whole triangle case
%        - K.CTP1,K.CTP2 for the part triangle case
%        (all are nCTs x nwt)
% Output: triplet on CTs: vv1 for ij, vv2 for k(i)k(j) (column-arrays)

iPs=CT.iPs; typeCTs=CT.type;nodeCTs=CT.nodes;
nCTs = size(CTs,2); % number of all not cut triangles in Omg1

% default K
if ~isempty(K)
    KW1=K.CTW1; KW2=K.CTW2; KP1=K.CTP1; KP2=K.CTP2; 
else
    dim=2; deg=pa.degN; 
    [wt,~] = getGaussQuad(dim,deg); 
    nwt = size(wt,2);
    Kdf = ones(nCTs,nwt);
    KW1=Kdf.CTW1; KW2=Kdf.CTW2; KP1=Kdf.CTP1; KP2=Kdf.CTP2;
end

ii = zeros(9*nCTs,1); jj = zeros(9*nCTs,1); 
vv1 = zeros(9*nCTs,1); vv2 = zeros(9*nCTs,1);
points = msh.p;


idx = 1;
for t=1:nCTs
    triangle = CTs(:,t); % info of triangle t-th
    iP1 = iPs(:,1,t); % 1st intersection point
    iP2 = iPs(:,2,t); % 2nd intersection point
    KKW1=KW1(t,:); KKW2=KW2(t,:); KKP1=KP1(t,:); KKP2=KP2(t,:);
    for i=1:3
        for j=1:3
            ii(idx) = CTs(i,t);
            jj(idx) = CTs(j,t);
            if typeCTs(t)==2 % 1 node in Omg2, 2 nodes in Omg1
                rV = points(:,nodeCTs.eachOmg2(1,t)); % the only node in Omg2
                vWhole1 = getTriplePPWhole(msh,pa,triangle,j,i,KKW1);
                vPart1 = getTriplePPPart(msh,pa,triangle,j,i,iP1,iP2,rV,KKP1);
                vv1(idx) = vWhole1 - vPart1;
                vv2(idx) = getTriplePPPart(msh,pa,triangle,j,i,iP1,iP2,rV,KKP2);
                idx = idx+1;
            else % typeCT = 0 or 4
                rV = points(:,nodeCTs.eachOmg1(1,t)); % the only node in Omg1
                vWhole2 = getTriplePPWhole(msh,pa,triangle,j,i,KKW2);
                vPart2 = getTriplePPPart(msh,pa,triangle,j,i,iP1,iP2,rV,KKP2);
                vv2(idx) = vWhole2 - vPart2;
                vv1(idx) = getTriplePPPart(msh,pa,triangle,j,i,iP1,iP2,rV,KKP1);
                idx = idx+1;
            end
        end
    end
end

end
