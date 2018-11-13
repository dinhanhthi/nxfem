function valG = getGij_PartTri(triangle,areaT,ii,jj,kk,iP1,iP2,remainVertex,dim,deg,msh)
% Find integrate of g(Phi_i,Phi_j) on a part of the single triangle (triangle-shape one)
% Input: - which phi_i (ii) & phi_j (jj)
%        - what is g(.,.)? (kk)
%        - 2 intersection points (iP1,iP2)
%        - only 1 remaining vertex: vInOmg1
%        - dimension dim and degree of quadrature
% Output: the value of integrate on each part of this triangle wrt nodes ii and jj 

points = msh.p;

[wt,pt] = getGaussQuad(dim,deg); % quadrature weights and points
nwt = size(wt,2); % number of Gaussian points

v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3

[xiP1h,yiP1h] = getCoorRef(iP1,v1,v2,v3); % cut point 1 in reference coordinate
[xiP2h,yiP2h] = getCoorRef(iP2,v1,v2,v3); % cut point 2 in reference coordinate
[xRvh,yRvh] = getCoorRef(remainVertex,v1,v2,v3); % remaining vertex in ref coor

% ref-triangle's part's area
areaTHp = getAreaTri([xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);


valG = 0;
for k=1:nwt
    % point in intermediate coordinate wrt the Gaussian point in reference coordiante
    [xHk,yHk] = getCoorSTD(pt(:,k),[xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);
    % shape functions at quadrature point in intermediate coordinate
    [shapeFunction,~,~] = getP1shapes(xHk,yHk);
    valG = valG + 2*areaT*areaTHp*wt(k)*defG(shapeFunction(ii),shapeFunction(jj),kk);
end

end