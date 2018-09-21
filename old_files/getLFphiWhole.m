function Fti = getLFphiWhole(triangle,i,msh,pa,defF,sub)
% int_Omg f*phi on NCTs (whole triangle)
% mother: getLFphiNCTs.m
% Input: - which node i?
%        - which triangle? (with points to get vertices)
%        - f in which domain? sub
% Output: load vector's value wrt this triangle at node ii

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
    [xk,yk] = getCoorSTD(pt(:,k),v1,v2,v3);
    [shFu,~,~] = getP1shapes(pt(1,k),pt(2,k)); % N_i at quadrature points
    Fti = Fti + areaT*wt(k)*defF(xk,yk,sub,pa)*shFu(i);
end

end