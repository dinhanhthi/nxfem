function gradPhi = getGradPhi(Ts,msh)
% Find gradient (Grad) of function phi_i of the triangles
% State: wrong!!!!
% Input: triangles we want to find Grad of Phi
% Output: matrix contains grad of all Phi_i of all triangles we want
%         matrix: 2 coordinates of grad x 3 vertices x nTris
% not that the index of last component in gradPhi (nTris) is the index of
%           triangle in tris, not in the global triangles!!!

points = msh.p; nTs = size(Ts,2);
gradPhi = zeros(2,3,nTs); % guarantee size of 2x3xnTs

tmp(:,1,:) = points(:,Ts(3,:)) - points(:,Ts(2,:)); % k-j
tmp(:,2,:) = points(:,Ts(1,:)) - points(:,Ts(3,:)); % i-k
tmp(:,3,:) = points(:,Ts(2,:)) - points(:,Ts(1,:)); % j-i

areaT = 0.5*abs(-tmp(1,3,:).*tmp(2,2,:) + tmp(1,2,:).*tmp(2,3,:));

for t=1:nTs
    % tmp = vector of (a,b)
    gradPhi(1,:,t) = -tmp(2,:,t); % xGrad = -b
    gradPhi(2,:,t) = tmp(1,:,t); % yGrad = a
end

for i=1:3
    gradPhi(1,i,:) = 1/(2*areaT).*gradPhi(1,i,:); % x
    gradPhi(2,i,:) = 1/(2*areaT).*gradPhi(2,i,:); % y
end

end