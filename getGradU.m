function gradU = getGradU(u,tris,msh)
% Get gradient of function u in standard FEM on triangles tris
% Status: checked with simple vectors
% Input: - u: number of standards nodes x 1\
%        - which triangles are considered? tris
%        - mesh msh
% Output: matrix: 2 coordinates x 3 vertices x number of triangles tris

ntris = size(tris,2); % number of triangles tris
tmp = zeros(2,3,ntris);

% gradient of all basis functions on all triangles tris
gradPhi = getGradPhi(tris,msh); % 2 coor x 3 vertices x nTs

for t=1:ntris
    tmp(:,1,t) = u(tris(1,t))*gradPhi(:,1,t);
    tmp(:,2,t) = u(tris(2,t))*gradPhi(:,2,t);
    tmp(:,3,t) = u(tris(3,t))*gradPhi(:,3,t);
end

gradU = tmp;

end