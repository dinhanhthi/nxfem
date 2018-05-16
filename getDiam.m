function hT = getDiam(msh)
% Get the diameter (longest side) of each triangle in the mesh
% State: checked with numSeg=2 and 3 + checked hTmax with freefem++
% Input: msh (points and triangles)
% Output: vector of diam of all triangles: 1 x number of triangles

triangles = msh.t; points = msh.p;

nTs = size(triangles,2); % number of triangles
ed = zeros(3,nTs); % contains all edges' size of each triangle
ed(1,:) = ((points(1,triangles(2,:))-points(1,triangles(3,:))).^2 ... %x
    + (points(2,triangles(2,:))-points(2,triangles(3,:))).^2).^(0.5); %y
ed(2,:) = ((points(1,triangles(3,:))-points(1,triangles(1,:))).^2 ... %x
    + (points(2,triangles(3,:))-points(2,triangles(1,:))).^2).^(0.5); %y
ed(3,:) = ((points(1,triangles(1,:))-points(1,triangles(2,:))).^2 ... %x
    + (points(2,triangles(1,:))-points(2,triangles(2,:))).^2).^(0.5); %y

hT = max(ed);
end