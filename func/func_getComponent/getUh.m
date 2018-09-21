function uh = getUh(wh,vh,kk1,kk2,msh)
% find uh = wh - kk*vh
% all uh, vh, wh are in NXFEM
% kk is different in each subdomain
% Input: wh, vh, kk
% Output: uh in NXFEM

newNodes = msh.newNodes;

uh = zeros(msh.ndof,1);
uh(msh.node.omg1) = wh(msh.node.omg1) - kk1*vh(msh.node.omg1);
uh(msh.node.omg2.all) = wh(msh.node.omg2.all); % v=0 on Omg2
uh(msh.node.onG) = wh(msh.node.onG); % v(onG)=0
uh(newNodes(msh.node.CT.all)) = wh(newNodes(msh.node.CT.all))...
                                - kk1*vh(newNodes(msh.node.CT.all));
    % note that there are nodes in msh.node.CT.all but in Omg2
    % thus we need the next line to find again uh at these nodes
uh(newNodes(msh.node.CT.omg2)) = wh(newNodes(msh.node.CT.omg2));
end