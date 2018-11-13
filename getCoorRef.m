function [xx,yy] = getCoorRef(pSTD,v1,v2,v3)
% Get corresponding point in the reference triangle (0-1-1)
% Input: - point in standard triangle
%        - 3 vertices v1 v3 v3 of the triangle
% Output: point in reference triangle coordinate

areaT = getAreaTri(v1,v2,v3); % area of triangle
xx = 1/(2*areaT)*((pSTD(1)-v1(1))*(v3(2)-v1(2))-(pSTD(2)-v1(2))*(v3(1)-v1(1)));
yy = 1/(2*areaT)*(-(pSTD(1)-v1(1))*(v2(2)-v1(2))+(pSTD(2)-v1(2))*(v2(1)-v1(1)));
end