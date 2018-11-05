function uNorVec = getUnitNV(pointA,pointB)
% Find the unit normal vector at a point on AB
% Input: Two endpoints A,B in that order with their coordinates (x,y)
% Note that, (A,B) is different from (B,A), the order is very important!
% Output: unit normal vector to that segment, this vector is always on the
%       left, that's why the order of A and B is very important

% the lenght of AB
lenAB = sqrt((pointB(1)-pointA(1))^2+(pointB(2)-pointA(2))^2);

% vector AB/|AB|
vec = (1/lenAB)*(pointB - pointA);

uNorVec(1) = -vec(2);
uNorVec(2) = vec(1);
end