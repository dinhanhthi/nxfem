function uStd = interNX2STD(uNX,msh)
% Interpolation from nxfem to std fem
% Input: solution in nxfem
% Output: solution in std fem (column array)

bNodesStd = unique([msh.e(1,:) msh.e(2,:)]); % standard boundary nodes
bNodesOnGam = intersect(msh.node.onG,bNodesStd); 
    % nodes both on Gam and on boundary

uStd(msh.node.omg1) = uNX(msh.node.omg1);
uStd(msh.node.omg2.notCT) = uNX(msh.node.omg2.notCT);
uStd(msh.node.CT.iomg2) = uNX(msh.newNodes(msh.node.CT.iomg2));
uStd(msh.node.CT.onG) = uNX(msh.newNodes(msh.node.CT.onG)) + uNX(msh.node.CT.onG);
uStd(bNodesOnGam) = uNX(bNodesOnGam);
uStd = uStd'; % transform to column-array for pdesurf

end