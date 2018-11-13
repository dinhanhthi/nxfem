function gamh = getGamh(iPs)
% Find matrix contains all coordinates of points of Gam_h (interface in Vh)
% Purpose: to plot together with the solution/mesh
% Input: intersection points
% Output: matrix 2 coordinates x number of cut points (no duplicate)

nCTs = size(iPs,3); % number of cut triangles
cPs = zeros(2,2*nCTs); % store the cut points
idx = 1;
for t=1:nCTs
    cPs(:,idx) = iPs(:,1,t);
    cPs(:,idx+1) = iPs(:,2,t);
    idx = idx+2;
end

% remove the duplicate points
cPs = unique(cPs','rows'); % nCPs x 2
cPs = cPs'; % convert to 2 x nCPs

gamh = cPs;
end