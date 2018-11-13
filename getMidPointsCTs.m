function mp = getMidPointsCTs(iPs)
% Find the mid points of gam_h on each cut triangles
% Purpose: plot the unit normal vector
% Input: intersection points
% Output: matrix 2 coordinates x nCTs

nCTs = size(iPs,3); % number of cut triangles
tmp = zeros(1,nCTs);

for t=1:nCTs
    tmp(1,t) = 0.5*(iPs(1,1,t)+iPs(1,2,t));
    tmp(2,t) = 0.5*(iPs(2,1,t)+iPs(2,2,t));
end

mp = tmp;
end