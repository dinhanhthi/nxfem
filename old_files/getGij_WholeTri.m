function valGij = getGij_WholeTri(triangle,ii,jj,pa,msh)
% Find integrate of Phi_i*Phi_j on the single whole triangle
% Input: - which phi_i (ii) & phi_j (jj)
%        - dimension dim and degree of quadrature
% Output: the value of integrate on this triangle wrt nodes ii and jj 

dim=2; deg=pa.degP2D; %Gaussian quadrature points in 2D
[wt,pt] = getGaussQuad(dim,deg); % quadrature weights and points
nwt = size(wt,2); % number of Gaussian points
points = msh.p;

v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3
areaT = getAreaTri(v1,v2,v3); % area of triangle  

valGij = 0;
for k=1:nwt
    [shFu,~,~] = getP1shapes(pt(1,k),pt(2,k)); % shape functions at quadrature points
    valGij = valGij + areaT*wt(k)*shFu(ii)*shFu(jj);
end

end