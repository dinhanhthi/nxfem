function [iN,bN] = getibNodes(msh)
% Get boundary and inner nodes for NXFEM
% Input:
% Output: inner nodes (iN), boundary nodes (bN) (all in row arrays)
%           (including the new nodes)
% DIFFERENT FROM getNodes or CT.nodes, NODES IN THIS FILE CONTAIN THE NEWNODES TOO!!!!

newNodes = msh.newNodes; edges = msh.e;
nodeCTs = msh.node.CT.all;

%% boundary nodes
bN.std = unique([edges(1,:) edges(2,:)]); % standard boundary nodes
bN.CT = intersect(bN.std,nodeCTs); % boundary nodes around interface
tmp = newNodes(bN.CT); % NEW boundary nodes around interface
bN.CT = [bN.CT,tmp']; % boundary nodes in CTs (new+std)
bN.all = [bN.std,tmp']; % all boundary nodes (std + new)


%% bN depends on ID of the edges
% note that, there are nodes at the corners (i.e, they are both on 2 edges)
tmp = find(edges(5,:)==1); % edge 1 (bottom)
bN.e1 = unique([edges(1,tmp) edges(2,tmp)]); % get std nodes
tmp = intersect(bN.e1,nodeCTs); % nodes in CTs
tmp = newNodes(tmp);
bN.e1 = [bN.e1,tmp']; % std+new

tmp = find(edges(5,:)==2); % edge 2 (right)
bN.e2 = unique([edges(1,tmp) edges(2,tmp)]); % get std nodes
tmp = intersect(bN.e2,nodeCTs); % nodes in CTs
tmp = newNodes(tmp);
bN.e2 = [bN.e2,tmp']; % std+new

tmp = find(edges(5,:)==3); % edge 3 (top)
bN.e3 = unique([edges(1,tmp) edges(2,tmp)]); % get std nodes
tmp = intersect(bN.e3,nodeCTs); % nodes in CTs
tmp = newNodes(tmp);
bN.e3 = [bN.e3,tmp']; % std+new

tmp = find(edges(5,:)==4); % edge 4 (left)
bN.e4 = unique([edges(1,tmp) edges(2,tmp)]); % get std nodes
tmp = intersect(bN.e4,nodeCTs); % nodes in CTs
tmp = newNodes(tmp);
bN.e4 = [bN.e4,tmp']; % std+new


%% inner nodes
iN.std = setdiff(1:msh.nStd,bN.std); % standard inner nodes
iN.CTs = intersect(iN.std,nodeCTs); % inner nodes around interface
tmp = newNodes(iN.CTs); % NEW inner nodes around interface
iN.CTs = [iN.CTs,tmp']; % inner nodes in CTs (new+std)
iN.all = [iN.std,tmp']; % all inner nodes (std + new)

end