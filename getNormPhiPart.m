function normPart = getNormPhiPart(triangle,areaT,ii,iP1,iP2,rV,dim,deg,msh)
% Find norm of phi_ii on a part of triangle (triangle-shape one)
% state: checked with mdl=4, tH=3, regu=0, nSeg=5
% This function get the idea from the function getLoadPartTri.m
% Input: - which node ii?
%        - which triangle? (with points to get vertices)
%        - triangle's area: areaT
%        - ref-triangle's part's area: areaTp
%        - 2 intersection points: iP1, iP2
%        - only 1 remaining vertex: rV
%        - dimension dim and degree of quadrature
% Output: - norm of phi_ii on the triangle part

points=msh.p;

[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points

v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3

[xiP1h,yiP1h] = getCoorRef(iP1,v1,v2,v3); % cut point 1 in reference coordinate
[xiP2h,yiP2h] = getCoorRef(iP2,v1,v2,v3); % cut point 2 in reference coordinate
[xRvh,yRvh] = getCoorRef(rV,v1,v2,v3); % remaining vertex in ref coor

% ref-triangle's part's area
areaTHp = getAreaTri([xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);

normPart = 0;
for k=1:nwt
    % point in intermediate coordinate wrt the Gaussian point in reference coordiante
    [xHk,yHk] = getCoorSTD(pt(:,k),[xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);
    % shape function N_i at quadrature point in intermediate coordinate
    [shapeFunction,~,~] = getP1shapes(xHk,yHk);
    normPart = normPart + 2*areaT*areaTHp*wt(k)*shapeFunction(ii)*shapeFunction(ii);
end

end