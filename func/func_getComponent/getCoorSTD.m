function [xx,yy] = getCoorSTD(pRef,v1,v2,v3)
% Get corresponding point in the standard triangle (v1-v2-v3)
% Input: - point in reference triangle (quadrature points)
%        - 3 vertices v1 v3 v3 of the triangle
% Output: point in standard triangle coordinate

xx = v1(1)*(1-pRef(1)-pRef(2))+v2(1)*pRef(1)+v3(1)*pRef(2);
yy = v1(2)*(1-pRef(1)-pRef(2))+v2(2)*pRef(1)+v3(2)*pRef(2);
end