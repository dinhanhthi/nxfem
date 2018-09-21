function [FCT1,FCT2] = getLoaduCTs(CTs,iPs,typeCTs,nodeCTs,msh,pa,defF,uold,wS)
% Get int_CTs (fu-wg(u))phi for equation of u
% State: 
% Input: - interPoints: 2 coor x 2 cut points x nCTs
%        - CTs:  5 x nCT (3 first line contain vertices)
%        - points: 2 x nP
%        - nodeCTs : nodes on each subdomain/on the interface of each triangle in CTs 
%        - wS, uold (in seperated subdomain)
% Output: load vector defined on cut triangles

% pa.degN: Gaussian quadrature points (for complicated function)
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
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
            Fwhole1 = getLoaduWhole(tri,1,i,dim,deg,msh,pa,defF,uoldT1,wST1);
            Fpart1 = getLoaduPart(tri,i,iP1,iP2,rV,dim,deg,1,msh,pa,defF,uoldT1,wST1);
            ff1(idx) = Fwhole1 - Fpart1;
            ff2(idx) = getLoaduPart(tri,i,iP1,iP2,rV,dim,deg,2,msh,pa,defF,uoldT2,wST2);
            idx = idx+1;
        else % typeCT = 0 or 4
            rV = points(:,nodeCTs.eachOmg1(1,t)); % the only vertex in Omg1
            ff1(idx) = getLoaduPart(tri,i,iP1,iP2,rV,dim,deg,1,msh,pa,defF,uoldT1,wST1);
            Fwhole2 = getLoaduWhole(tri,2,i,dim,deg,msh,pa,defF,uoldT2,wST2);
            Fpart2 = getLoaduPart(tri,i,iP1,iP2,rV,dim,deg,2,msh,pa,defF,uoldT2,wST2);
            ff2(idx) = Fwhole2 - Fpart2;
            idx = idx+1;
        end % end if typeCT
    end % end for vertices
end % end for nCT

FCT1 = accumarray(ii,ff1);
FCT2 = accumarray(ii,ff2);

end