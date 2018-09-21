function intWhole = getTriplegPPWhole(triangle,uoldT,i,j,dim,deg,...
                        msh,typeG,kk)
% Find int_triangle(kk*g(uold)*phi*phi)on whole triangle
% state: 
% Input: - which node i,j
%        - triangle's area
%        - coefficient kk (any)
%        - uoldT: value of uold at vertices of triangle T
%        - dimension dim and degree of quadrature
% Output: - the value on whole triangle

[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points
points = msh.p;

v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3
areaT = getAreaTri(v1,v2,v3); % area of triangle  

intWhole=0;
for k=1:nwt
    [shFu,~,~] = getP1shapes(pt(1,k),pt(2,k)); % N_i at quadrature points
    vuold = uoldT(1)*shFu(1)+uoldT(2)*shFu(2)+uoldT(3)*shFu(3);
    guold = defG(vuold,typeG); % g(uold)
    intWhole = intWhole + kk*areaT*wt(k)*guold*shFu(i)*shFu(j);
end

end