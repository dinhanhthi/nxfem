function unx = interSTD2NX(ustd,msh)
% Interpolation from std fem to nxfem
% Input: solution in standard fem
% Output: solution in nxfem (column array)

unx = zeros(msh.ndof,1); % column-array
unx(msh.node.std) = ustd(msh.node.std);
unx(msh.newNodes(msh.node.CT.all)) = ustd(msh.node.CT.all);

end