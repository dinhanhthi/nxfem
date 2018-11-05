function nodes = getNodes(tris,nodeCTs,msh,phi,pa)
% get all necessary nodes of the mesh
% Input:
% Output: structure var contains all nodes (std nodes)


%% NORMAL NODES
nodes.CT.onG = nodeCTs.onG; % nodes of CTs on Gam
nodes.CT.iomg2 = nodeCTs.Omg2; % nodes in Omg2 (not on interface) of CTs
nodes.CT.all = nodeCTs.all; % nodes in CTs, row-array
nodes.CT.omg2 = [nodeCTs.Omg2,nodeCTs.onG]; % nodes in Omg2 and on G of CTs, row-array

nodes.std = 1:msh.nStd; % standard nodes
nodes.iomg1 = nodes.std((phi<0)&(abs(phi)>pa.tol)); % all std nodes inside Omg1
nodes.onG = nodes.std(abs(phi)<pa.tol); % all std nodes on Gam
nodes.iomg2 = nodes.std((phi>0)&(abs(phi)>pa.tol)); % all std nodes inside Omg2
    
nodes.omg1 = unique(tris.NCTs1(1:3,:));
nodes.omg1 = union(nodes.omg1,nodeCTs.onG);
    % nodes inside Omg1 and on Gam, column-array
    % there still nodes on Gam in CTs which not belongs to 
    %   NCTs, that's why we need 2nd line
nodes.omg2.all = unique(tris.NCTs2(1:3,:)); 
    % column-array, there still nodes on Gam in CTs not belongs to NCTs
nodes.omg2.all = union(nodes.omg2.all,nodeCTs.onG);
    % nodes inside Omg2 and on Gam, column-array
nodes.omg2.notCT = setdiff(nodes.omg2.all,nodeCTs.all);
    % nodes in Omg2 but not in CTs (colum-array)
    % there may be still some nodes on the interface

end