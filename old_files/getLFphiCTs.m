function [ii,ff1,ff2] = getLFphiCTs(CTs,iPs,typeCTs,nodeCTs,msh,pa,defF)
% int_Omg f*phi on CTs (load vector)
% State: alike getLoadCTs, diff at 16th digit after ","
% Input: - interPoints: 2 coor x 2 cut points x nCTs
%        - CTs:  5 x nCT (3 first line contain vertices)
%        - points: 2 x nP
%        - nodeCTs : nodes on each subdomain/on the interface of each triangle in CTs 
% Output: ii (nodes); ff1, ff2 (values at nodes). column arrays

points=msh.p;

nCTs = size(CTs,2); % number of cut triangles
ii = zeros(3*nCTs,1); ff1 = zeros(3*nCTs,1); ff2 = zeros(3*nCTs,1);

idx=1;
for t=1:nCTs
    iP1 = iPs(:,1,t); % 1st intersection point
    iP2 = iPs(:,2,t); % 2nd intersection point
    tri = CTs(:,t);
    for i=1:3 % 3 vertices
        ii(idx) = CTs(i,t);
        if typeCTs(t)==2 % 1 node in Omg2, 2 nodes in Omg1
            rV = points(:,nodeCTs.eachOmg2(1,t)); % the only vertex in Omg2
            Fwhole1 = getLFphiWhole(tri,i,msh,pa,defF,1);
            Fpart1 = getLFphiPart(tri,i,iP1,iP2,rV,msh,pa,defF,1);
            ff1(idx) = Fwhole1 - Fpart1;
            ff2(idx) = getLFphiPart(tri,i,iP1,iP2,rV,msh,pa,defF,2);
            idx = idx+1;
        else % typeCT = 0 or 4
            rV = points(:,nodeCTs.eachOmg1(1,t)); % the only vertex in Omg1
            Fwhole2 = getLFphiWhole(tri,i,msh,pa,defF,2);
            Fpart2 = getLFphiPart(tri,i,iP1,iP2,rV,msh,pa,defF,2);
            ff2(idx) = Fwhole2 - Fpart2;
            ff1(idx) = getLFphiPart(tri,i,iP1,iP2,rV,msh,pa,defF,1);
            idx = idx+1;
        end % end if typeCT
    end % end for vertices
end % end for nCT

end