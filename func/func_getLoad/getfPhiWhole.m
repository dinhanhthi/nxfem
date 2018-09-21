function val = getfPhiWhole(msh,pa,triangle,vertex,P)
% int_Omg P*phi on NCTs/CTs (whole triangle)
% This file is used for general P (the size of P depends on number of
%       Gaussian points (pa.degN)
% Note: rewriten after 1st rewrite with inputParser (worst performance)
% State: checked with getLf
% Related: getfPhiNCTs, getfPhiCTs, getfPhiPart
% Input: - which vertex (1,2 or 3)
%        - which triangle? (with points to get vertices)
%        - component P: 1 x nwt
% Output: scalar value
    
points=msh.p;
% pa.degN: Gaussian quadrature points (for complicated function)
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points
v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3
areaT = getAreaTri(v1,v2,v3); % area of triangle

% default P
if isempty(P)
    P = ones(1,nwt);
end

val=0;
for k=1:nwt
    [shFu,~,~] = getP1shapes(pt(1,k),pt(2,k)); % N_i at quadrature points
    val = val + P(k)*areaT*wt(k)*shFu(vertex);
end

end
