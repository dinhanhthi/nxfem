function pGh = getPointsGamh(msh,pa,iPs)
% Get all points of Gamma_h (including intersection points and points
%   coinciding to the mesh's nodes.
% Note: it's different from CT.iPs because it includes also the nodes in
%   the fitted-mesh case.
% First used for plotting interface.
%   Note: DON'T USE because the order of points is arbitary while we need to plot the interface in some ordered points.
% Input: - iPs: intersection points: 2 coors x 2 cut points x nCTs
%        - msh: mesh data (including msh.node)
% Output: a vector pGh containing all points (2 coor x number of points)

A1 = []; A2 = [];

% collect all interface nodes lying on cut triangles
%   (there may be similarities)
if ~isempty(iPs)
    nCTs = size(iPs,3); % number of CTs
    A1 = reshape(iPs, 2, nCTs*2); % 2 x #nodes
end

% collect all nodes which interface goes through
if ~isempty(msh.node.onG)
    A2 = msh.p(:,msh.node.onG); % 2 x #nodes
end

pGh = [A1, A2];
pGh(abs(pGh)<pa.tol) = 0; % make all "very small" values be 0
pGh = unique(pGh','rows');
pGh = pGh';


end