function areaChildCT = getAreaChild(cutTriangles,areaCT,interPoints,...
                            idxOmg1CT,idxOmg2CT,typeCT,msh)
% Find area of a child of a cut triangle
% Input: cut triangles, phi, points and case of Omg1 or Omg2
% Output: 2 x nCT whose 1st row describes Omg1, 2nd row desribes Omg2

nCT = size(cutTriangles,2);
points = msh.p;

tmp = zeros(2,nCT); 
for i=1:nCT
    if typeCT(i)==2 % 2 nodes in Omg1, 1 node in Omg2
        tmp(2,i) = getAreaTri(interPoints(:,1,i),interPoints(:,2,i),points(:,idxOmg2CT(1,i)));
        tmp(1,i) = areaCT(i) - tmp(2,i);
    elseif typeCT(i)==4 % 2 nodes in Omg2, 1 node in Omg1
        tmp(1,i) = getAreaTri(interPoints(:,1,i),interPoints(:,2,i),points(:,idxOmg1CT(1,i)));
        tmp(2,i) = areaCT(i) - tmp(1,i);
    else % type==0, 1 node lies on the interface
        tmp(1,i) = getAreaTri(interPoints(:,1,i),interPoints(:,2,i),points(:,idxOmg1CT(1,i)));
        tmp(2,i) = areaCT(i) - tmp(1,i);
    end
end

areaChildCT = tmp;
end