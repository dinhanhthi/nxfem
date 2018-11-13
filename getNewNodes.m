function newNodes = getNewNodes(nodesAroundGamCTs,nNodes)
% Map current node to the newNode depending on the number of Node nNodes
% Input: nodes around the interface and total nodes of the mesh
% Output: a vector contain new numbering-nodes whose indices are nodes
%           around the interface

newNodes=NaN(1,max(nodesAroundGamCTs));
nNodesAroundGam = size(nodesAroundGamCTs,2);
newNodes(nodesAroundGamCTs)=nNodes+1:nNodes+nNodesAroundGam;
newNodes = newNodes'; % transform to column-array
end