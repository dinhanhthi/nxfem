function gP = getGrad(i,t,msh)
% Get grad of phi at vertex i of triangle t
% Status: checked!
% Input: - vertex i (from 1 to 3 of triangle t)
%        - triangle t trong msh.t
%        - mesh msh
% Output: vector 2 coordinates x 1

points = msh.p; triangles = msh.t;
% j = mod(i+1,3)+3*floor((3-mod(i+1,3))/3); % vertex j
% k = mod(i+2,3)+3*floor((3-mod(i+2,3))/3); % vertex k
j = -3/2*i^2+11/2*i-2; % vertex j
k = 3/2*i^2-13/2*i+8; % vertex k
v(:,1) = points(:,triangles(1,t)); % vertex 1 of triangle t
v(:,2) = points(:,triangles(2,t)); % vertex 2 of triangle t
v(:,3) = points(:,triangles(3,t)); % vertex 3 of triangle t
at = getAreaTri(v(:,1),v(:,2),v(:,3)); % area of triangle t

gP = zeros(2,1);
gP(1) = 1/(2*at)*(points(2,triangles(j,t))-points(2,triangles(k,t))); % grad_x
gP(2) = 1/(2*at)*(points(1,triangles(k,t))-points(1,triangles(j,t))); % grad_y
end