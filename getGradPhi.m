function gradPhi = getGradPhi(tris,msh)
% Find gradient (Grad) of function phi_i of the triangles
% State: checked with numSeg=2
% Input: triangles we want to find Grad of Phi
% Output: matrix contains grad of all Phi_i of all triangles we want
%         matrix: 2 coordinates of grad x 3 vertices x nTris
% not that the index of last component in gradPhi (nTris) is the index of
%           triangle in tris, not in the global triangles!!!

points = msh.p;

tmp(:,1,:) = points(:,tris(3,:)) - points(:,tris(2,:)); % k-j
tmp(:,2,:) = points(:,tris(1,:)) - points(:,tris(3,:)); % i-k
tmp(:,3,:) = points(:,tris(2,:)) - points(:,tris(1,:)); % j-i

areaT = 0.5*abs(-tmp(1,3,:).*tmp(2,2,:) + tmp(1,2,:).*tmp(2,3,:));

% tmp = vector of (a,b)
gradPhi(1,:,:) = -tmp(2,:,:); % xGrad = -b
gradPhi(2,:,:) = tmp(1,:,:); % yGrad = a

for i=1:3
    gradPhi(1,i,:) = 1/(2*areaT).*gradPhi(1,i,:); % x
    gradPhi(2,i,:) = 1/(2*areaT).*gradPhi(2,i,:); % y
end

end