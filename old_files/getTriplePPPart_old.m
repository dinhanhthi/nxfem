function intPart = getTriplePPPart(triangle,i,j,iP1,iP2,rV,pa,msh,K)
% Find int_NCTs(K*phi*phi) on part triangle (triangle shape one)
% This file is used for general K (the size of K depends on number of
%       Gaussian points (dim,deg)
% state: 
% Input: - which triangle
%        - which node i,j
%        - 2 intersection points: iP1, iP2
%        - remaining vertex: rV
%        - dimension dim and degree of quadrature
%        - component K
% Output: the value of integrate on each part of triangle wrt nodes i,j

points = msh.p;
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D (polynomial func)
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
    intPart = intPart + K(k)*2*areaT*areaTHp*wt(k)*shFu(i)*shFu(j);
end

end