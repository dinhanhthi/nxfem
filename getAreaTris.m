function areaTris = getAreaTris(tris,msh)
% Function to find area of many triangles at the same time
% Input: triangles
% Output: area of these triangles (1 x nTris)

points=msh.p;

% normal vectors
% ve is a 2 coordinates x nTriangles x 3 edges
ve(:,:,3) = points(:,tris(2,:)) - points(:,tris(1,:));
ve(:,:,1) = points(:,tris(3,:)) - points(:,tris(2,:));
ve(:,:,2) = points(:,tris(1,:)) - points(:,tris(3,:));

% area of the triangles
areaTris = 0.5*abs(ve(1,:,3).*ve(2,:,2) - ve(1,:,2).*ve(2,:,3)); 

end