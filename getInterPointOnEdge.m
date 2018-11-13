function interPointOnEdge = getInterPointOnEdge(point1,point2,phi1,phi2)
% Find the intersection point between phi and edge point1-pont2
% Input: - phi at 2 endpoint: 1x1
%        - 2 endpoints of the edge 2 coor x 1
% Output: intersection point

interPointOnEdge = 0.5*point1*(1-(phi1+phi2)/(phi1-phi2))... 
    + 0.5*point2*(1+(phi1+phi2)/(phi1-phi2));
end