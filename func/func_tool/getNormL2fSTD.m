function val = getNormL2fSTD(msh,pa,f)
% Get ||f(x,y)||_L2(\Omg) where f is a function handle
% First used in finding find norm of exact phi (level set) in Omg
% Input: - f(x,y,pa) : function handle
% Output: scalar value of norm
% ---------------------------------------------
% Important (3/11): this f in this case is NOT for the exact solution
%   uex(x,y,pa,sub), for this type of function (handle), cf.
%   getNormL2uexSTD.m
% ---------------------------------------------

triangles = msh.t; points = msh.p;
nTs = size(triangles,2);

dim=2; deg=pa.degN;
[wt,pt] = getGaussQuad(dim,deg);
nwt = size(wt,2);

val = 0;
for t = 1:nTs
    triangle = triangles(:,t);
    v1 = points(:,triangle(1)); % vertex 1
    v2 = points(:,triangle(2)); % vertex 2
    v3 = points(:,triangle(3)); % vertex 3
    areaT = getAreaTri(v1,v2,v3); % area of triangle
    tmp = 0;
    for k = 1:nwt
        [xk,yk] = getCoorSTD(pt(:,k),v1,v2,v3);
        tmp = tmp + wt(k)*f(xk,yk,pa)*f(xk,yk,pa);
    end
    
    val = val + areaT*tmp;
end

val = sqrt(val);

end