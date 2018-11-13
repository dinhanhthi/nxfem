function areaTri = getAreaTri(p1,p2,p3)
% Find the area of an arbitrary triangle
% Input: 3 vertices (x,y)
% Output: the area

% collect coordinates
xx = [p1(1), p2(1), p3(1)]; % x values
yy = [p1(2), p2(2), p3(2)]; % y values

% compute area
areaTri = polyarea(xx,yy);
end