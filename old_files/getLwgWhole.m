function Fti = getLwgWhole(triangle,i,pa,msh,uoldT,wST,typeG)
% Get int_NCTs wg(u)*phi on the whole triangle (for eqn u)
% Input: - which node ii?
%        - which triangle? (with points to get vertices)
%        - dimension dim and degree of quadrature
%        - uoldT and wST
% Output: load vector's value wrt this triangle at node i

points=msh.p;

% pa.degN: Gaussian quadrature points (for complicated function)
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points
v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3
areaT = getAreaTri(v1,v2,v3); % area of triangle

Fti=0;
for k=1:nwt
    [shFu,~,~] = getP1shapes(pt(1,k),pt(2,k)); % N_i at quadrature points
    vuold = uoldT(1)*shFu(1)+uoldT(2)*shFu(2)+uoldT(3)*shFu(3);
    guold = defG(vuold,typeG); % g(uold)
    vwS = wST(1)*shFu(1)+wST(2)*shFu(2)+wST(3)*shFu(3); % w
    Fti = Fti + vwS*guold*areaT*wt(k)*shFu(i);
end

end