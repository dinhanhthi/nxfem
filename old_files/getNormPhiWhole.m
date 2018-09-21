function normWhole = getNormPhiWhole(triangle,ii,pa,msh)
% Find norm of phi_ii on the whole triangle
% state: checked with mdl=4, tH=3, regu=0, nSeg=5
% Input: - which node ii?
%        - triangle's area
%        - dimension dim and degree of quadrature
% Output: - norm of phi_ii on the whole triangle

points = msh.p;
dim=2; deg=pa.degP2D; % Gaussian quadrature points in 2D
[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points


v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3
areaT = getAreaTri(v1,v2,v3);

normWhole=0;
for k=1:nwt
    [shFu,~,~] = getP1shapes(pt(1,k),pt(2,k)); % N_i at quadrature points
    normWhole = normWhole + areaT*wt(k)*shFu(ii)*shFu(ii);
end

end