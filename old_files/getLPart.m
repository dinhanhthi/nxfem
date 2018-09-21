function Fp = getLPart(triangle,i,iP1,iP2,rV,msh,pa,P)
% int_Omg P*phi on CTs (part triangle)
% This file is used for general P (the size of P depends on number of
%       Gaussian points (pa.degN)
% Input: - which node i?
%        - which triangle? (with points to get vertices)
%        - 2 intersection points: iP1, iP2
%        - only 1 remaining vertex: rV
%        - component P
% Output: load vector's value on a part of triangle at node i

points=msh.p;

% pa.degN: Gaussian quadrature points (for complicated function)
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points

v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3
areaT = getAreaTri(v1,v2,v3); % area of triangle

[xiP1h,yiP1h] = getCoorRef(iP1,v1,v2,v3); % cut point 1 in ref coor
[xiP2h,yiP2h] = getCoorRef(iP2,v1,v2,v3); % cut point 2 in ref coor
[xRvh,yRvh] = getCoorRef(rV,v1,v2,v3); % remaining vertex in ref coor
% ref-triangle's part's area
areaTHp = getAreaTri([xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);

Fp = 0;
for k=1:nwt
    % point in intermediate coordinate wrt the Gaussian point in ref coor
    [xHk,yHk] = getCoorSTD(pt(:,k),[xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);
    % shape function N_i at quadrature point in intermediate coordinate
    [shFu,~,~] = getP1shapes(xHk,yHk);
    Fp = Fp + P(k)*2*areaT*areaTHp*wt(k)*shFu(i);
end

end