function areaChildCT = getAreaChild(CTs,iPs,idxOmg1CT,idxOmg2CT,typeCT,msh)
% Find area of a child of a cut triangle
% Input: cut triangles, phi, points and case of Omg1 or Omg2
% Output: 2 x nCT whose 1st row describes Omg1, 2nd row desribes Omg2

nCT = size(CTs,2);
points = msh.p;

tmp = zeros(2,nCT); 
for t=1:nCT
    areaCT = getAreaTri(points(:,CTs(1,t)),points(:,CTs(2,t)),points(:,CTs(3,t)));
    if typeCT(t)==2 % 2 nodes in Omg1, 1 node in Omg2
        tmp(2,t) = getAreaTri(iPs(:,1,t),iPs(:,2,t),points(:,idxOmg2CT(1,t)));
        tmp(1,t) = areaCT - tmp(2,t);
    elseif typeCT(t)==4 % 2 nodes in Omg2, 1 node in Omg1
        tmp(1,t) = getAreaTri(iPs(:,1,t),iPs(:,2,t),points(:,idxOmg1CT(1,t)));
        tmp(2,t) = areaCT - tmp(1,t);
    else % type==0, 1 node lies on the interface
        tmp(1,t) = getAreaTri(iPs(:,1,t),iPs(:,2,t),points(:,idxOmg1CT(1,t)));
        tmp(2,t) = areaCT - tmp(1,t);
    end
end

areaChildCT = tmp;
end