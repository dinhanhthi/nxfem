function val = getNormL2GfSTD(msh,pa,fx,fy)
% Get ||grad f(x,y)||_L2 where f is a function handle
% IMPORTANT: This is not a general cases, need to be coded later!!
% First used in main_nxfem to find ENorm
% Input: - fx: \partial_x of f, function handle (x,y,pa)
%        - fy: \partial_y of f, function handle (x,y,pa)

nfx = getNormL2fSTD(msh,pa,fx);
nfy = getNormL2fSTD(msh,pa,fy);

val = nfx^2 + nfy^2;
val = sqrt(val);

end