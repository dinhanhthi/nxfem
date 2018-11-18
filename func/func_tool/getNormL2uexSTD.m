function val = getNormL2uexSTD(msh,pa,uex,defPhi)
% Get ||uex(x,y)||_L2(\Omg) where uex is a function handle
% First used in finding find norm of uex in Omg
% Input: - uex(x,y,pa,sub) : function handle
%        - defPhi(x,y,pa): function handle
% Output: scalar value of norm
% ---------------------------------------------
% Important (3/11): this file is diff (but not diff) from getNormL2STD, the
%   only diff is on the parser parameters uex(x,y,pa,sub) whereas f(x,y,pa).
%   In function uex, there is uex1, uex2 and we don't know how to determine
%   xk,yk in which subdomain?
% ---------------------------------------------

% IDEA (for future code): from values of phi, we can determind the
%   intersections. We can determine the position of xk,yk based on these
%   intersection.

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
        if abs(defPhi(xk,yk,pa))<pa.tol % in Omg1
            sub=1;
        else % in Omg2
            sub=2;
        end
        tmp = tmp + wt(k)*uex(xk,yk,pa,sub)*uex(xk,yk,pa,sub);
    end
    
    val = val + areaT*tmp;
end

val = sqrt(val);

end