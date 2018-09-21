function intWhole = getTriplePPWhole(triangle,i,j,pa,msh,K)
% Find int_NCTs(K*phi*phi) on whole triangle
% This file is used for general K (the size of K depends on number of
%       Gaussian points (dim,deg)
% state: 
% Input: - which node i,j
%        - dimension dim and degree of quadrature
%        - component K: nwt x 1
% Output: - the value on whole triangle

dim=2; deg=pa.degN; % Gaussian quadrature points in 2D (polynomial func)
[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points
points = msh.p;

v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3
areaT = getAreaTri(v1,v2,v3); % area of triangle 

if isempty(K)
    K = ones(nwt,1);
end

intWhole=0;
for k=1:nwt
    [shFu,~,~] = getP1shapes(pt(1,k),pt(2,k)); % N_i at quadrature points
    intWhole = intWhole + K(k)*areaT*wt(k)*shFu(i)*shFu(j);
end

end