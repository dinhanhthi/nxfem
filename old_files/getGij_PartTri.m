function valG = getGij_PartTri(triangle,ii,jj,iP1,iP2,rV,pa,msh)
% Find integrate of Phi_i*Phi_j on a part of the single triangle (triangle-shape one)
% Input: - which phi_i (ii) & phi_j (jj)
%        - 2 intersection points (iP1,iP2)
%        - only 1 remaining vertex: vInOmg1
%        - dimension dim and degree of quadrature
% Output: the value of integrate on each part of this triangle wrt nodes ii and jj 

points = msh.p;

dim=2; deg=pa.degP2D; %Gaussian quadrature points in 2D
[wt,pt] = getGaussQuad(dim,deg); % quadrature weights and points
nwt = size(wt,2); % number of Gaussian points

v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3
areaT = getAreaTri(v1,v2,v3);

[xiP1h,yiP1h] = getCoorRef(iP1,v1,v2,v3); % cut point 1 in reference coordinate
[xiP2h,yiP2h] = getCoorRef(iP2,v1,v2,v3); % cut point 2 in reference coordinate
[xRvh,yRvh] = getCoorRef(rV,v1,v2,v3); % remaining vertex in ref coor

% ref-triangle's part's area
areaTHp = getAreaTri([xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);


valG = 0;
for k=1:nwt
    % point in intermediate coordinate wrt the Gaussian point in ref coor
    [xHk,yHk] = getCoorSTD(pt(:,k),[xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);
    % shape functions at quadrature point in intermediate coordinate
    [shFu,~,~] = getP1shapes(xHk,yHk);
    valG = valG + 2*areaT*areaTHp*wt(k)*shFu(ii)*shFu(jj);
end

end