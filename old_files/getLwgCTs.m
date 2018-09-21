function [ii,ff1,ff2] = getLwgCTs(CTs,iPs,typeCTs,nodeCTs,msh,pa,uold,wS,typeG)
% Get int_CTs (fu-wg(u))phi for equation of u
% State: 
% Input: - interPoints: 2 coor x 2 cut points x nCTs
%        - CTs:  5 x nCT (3 first line contain vertices)
%        - points: 2 x nP
%        - nodeCTs : nodes on each subdomain/on the interface of each triangle in CTs 
%        - wS, uold (in seperated subdomain)
% Output: ii (nodes); ff1, ff2 (values at nodes). column arrays

points=msh.p;

nCTs = size(CTs,2); % number of cut triangles
ii = zeros(3*nCTs,1); ff1 = zeros(3*nCTs,1); ff2 = zeros(3*nCTs,1);

idx=1;
for t=1:nCTs
    uoldT1(1:3) = uold.ct1(CTs(1:3,t)); % uold's values at vertices of triangle
    wST1(1:3) = wS.ct1(CTs(1:3,t));
    uoldT2(1:3) = uold.ct2(CTs(1:3,t)); % uold's values at vertices of triangle
    wST2(1:3) = wS.ct2(CTs(1:3,t));
    iP1 = iPs(:,1,t); % 1st intersection point
    iP2 = iPs(:,2,t); % 2nd intersection point
    tri = CTs(:,t);
    for i=1:3 % 3 vertices
        ii(idx) = CTs(i,t);
        if typeCTs(t)==2 % 1 node in Omg2, 2 nodes in Omg1
            rV = points(:,nodeCTs.eachOmg2(1,t)); % the only vertex in Omg2
            Fwhole1 = getLwgWhole(tri,i,pa,msh,uoldT1,wST1,typeG);
            Fpart1 = getLwgPart(tri,i,iP1,iP2,rV,pa,msh,uoldT1,wST1,typeG);
            ff1(idx) = Fwhole1 - Fpart1;
            ff2(idx) = getLwgPart(tri,i,iP1,iP2,rV,pa,msh,uoldT2,wST2,typeG);
            idx = idx+1;
        else % typeCT = 0 or 4
            rV = points(:,nodeCTs.eachOmg1(1,t)); % the only vertex in Omg1
            ff1(idx) = getLwgPart(tri,i,iP1,iP2,rV,pa,msh,uoldT1,wST1,typeG);
            Fwhole2 = getLwgWhole(tri,i,pa,msh,uoldT2,wST2,typeG);
            Fpart2 = getLwgPart(tri,i,iP1,iP2,rV,pa,msh,uoldT2,wST2,typeG);
            ff2(idx) = Fwhole2 - Fpart2;
            idx = idx+1;
        end % end if typeCT
    end % end for vertices
end % end for nCT

end