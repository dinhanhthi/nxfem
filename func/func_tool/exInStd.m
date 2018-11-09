function ex = exInStd(defEx,msh,pa)
% exact solution in std fem
% Input: defition of exact solution via function defEx
% Output: ustd in std fem

points = msh.p;
ex = zeros(msh.nStd,1); % column-array
ex(msh.node.iomg1) = defEx(points(1,msh.node.iomg1),...
                        points(2,msh.node.iomg1),1,pa); % nodes inside Omg1
ex(msh.node.onG) = defEx(points(1,msh.node.onG),...
                        points(2,msh.node.onG),1,pa); % nodes on Gam
ex(msh.node.iomg2) = defEx(points(1,msh.node.iomg2),...
                        points(2,msh.node.iomg2),2,pa); % nodes inside Omg2

end