function [ii,ff1,ff2] = getLCTs(CTs,CT,msh,pa,P)
% int_Omg P*phi on CTs (load vector)
% This file is used for general P (the size of P depends on number of
%       Gaussian points (pa.degN)
% State: alike getLoadCTs, diff at 16th digit after ","
% Input: - interPoints: 2 coor x 2 cut points x nCTs
%        - CTs:  5 x nCT (3 first line contain vertices)
%        - points: 2 x nP
%        - nodeCTs : nodes on each subdomain/on the interface of each triangle in CTs
%        - P1,P2 for the whole/part case
% Output: ii (nodes); ff1, ff2 (values at nodes). column arrays

iPs=CT.iPs; typeCTs=CT.type; nodeCTs=CT.nodes;
points=msh.p;
PW1=P.CTW1; PW2=P.CTW2; PP1=P.CTP1; PP2=P.CTP2;

nCTs = size(CTs,2); % number of cut triangles
ii = zeros(3*nCTs,1); ff1 = zeros(3*nCTs,1); ff2 = zeros(3*nCTs,1);

idx=1;
for t=1:nCTs
    iP1 = iPs(:,1,t); % 1st intersection point
    iP2 = iPs(:,2,t); % 2nd intersection point
    tri = CTs(:,t);
    PPW1=PW1(t,:); PPW2=PW2(t,:); PPP1=PP1(t,:); PPP2=PP2(t,:);
    for i=1:3 % 3 vertices
        ii(idx) = CTs(i,t);
        if typeCTs(t)==2 % 1 node in Omg2, 2 nodes in Omg1
            rV = points(:,nodeCTs.eachOmg2(1,t)); % the only vertex in Omg2
            Fwhole1 = getLWhole(tri,i,msh,pa,PPW1);
            Fpart1 = getLPart(tri,i,iP1,iP2,rV,msh,pa,PPP1);
            ff1(idx) = Fwhole1 - Fpart1;
            ff2(idx) = getLPart(tri,i,iP1,iP2,rV,msh,pa,PPP2);
            idx = idx+1;
        else % typeCT = 0 or 4
            rV = points(:,nodeCTs.eachOmg1(1,t)); % the only vertex in Omg1
            Fwhole2 = getLWhole(tri,i,msh,pa,PPW2);
            Fpart2 = getLPart(tri,i,iP1,iP2,rV,msh,pa,PPP2);
            ff2(idx) = Fwhole2 - Fpart2;
            ff1(idx) = getLPart(tri,i,iP1,iP2,rV,msh,pa,PPP1);
            idx = idx+1;
        end % end if typeCT
    end % end for vertices
end % end for nCT

end