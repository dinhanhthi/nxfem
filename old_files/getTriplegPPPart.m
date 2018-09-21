function intPart = getTriplegPPPart(triangle,uoldT,i,j,iP1,iP2,rV,...
                                            dim,deg,typeG,msh,kk)
% Find int_triangle(kk*g(uold)*phi*phi)on a part of triangle (triangle
%               shape one)
% state: 
% Input: - which node i,j
%        - 2 intersection points: iP1, iP2
%        - remaining vertex: rV
%        - info of triangle: triangle
%        - coefficient kk (any)
%        - uoldT: value of uold at vertices of triangle T
%        - dimension dim and degree of quadrature
% Output: the value of integrate on each part of triangle wrt nodes i,j

points = msh.p;
[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points

v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3
areaT = getAreaTri(v1,v2,v3); % area of triangle  

[xiP1h,yiP1h] = getCoorRef(iP1,v1,v2,v3); % iP1 in reference coordinate
[xiP2h,yiP2h] = getCoorRef(iP2,v1,v2,v3); % iP2 in reference coordinate
[xRvh,yRvh] = getCoorRef(rV,v1,v2,v3); % remaining vertex in ref coor

% ref-triangle's part's area
areaTHp = getAreaTri([xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);

intPart = 0;
for k=1:nwt
    % point in intermediate coordinate wrt the Gaussian point in ref coor
    [xHk,yHk] = getCoorSTD(pt(:,k),[xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);
    % shape functions at quadrature point in intermediate coordinate
    [shFu,~,~] = getP1shapes(xHk,yHk);
    vuold = uoldT(1)*shFu(1)+uoldT(2)*shFu(2)+uoldT(3)*shFu(3);
    guold = defG(vuold,typeG); % g(uold)
    intPart = intPart + 2*areaT*areaTHp*kk*wt(k)*guold*shFu(i)*shFu(j);
end

end