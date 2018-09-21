function val = getfWhole(msh,pa,triangle,P)
% int_T P dx on whole triangle
% Note: the idea is the same with getfPhiWhole but there is no "phi"
% Status: 
% Related: - first used in getTripleGGVelo.m
%          - borrowed idea from getfPhiWhole.m
% Input: - which triangle? (with points to get vertices)
%        - component P: 1 x nwt
% Output: scalar value


points = msh.p; % points of the mesh
% pa.degN: Gaussian quadrature points (for complicated function)
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points

% default P
if isempty(P)
    P = ones(1,nwt);
end

v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3
areaT = getAreaTri(v1,v2,v3); % area of triangle


%% =======================================================================
% Get fPhiWhole
%=========================================================================
val=0;
for k=1:nwt
    [shFu,~,~] = getP1shapes(pt(1,k),pt(2,k)); % N_i at quadrature points
    [xk,yk] = getCoorSTD(pt(:,k),v1,v2,v3);
    val = val + rP.fw(vwS)*rP.gu(vuS)*rP.h(xk,yk,pa)*areaT*wt(k);
end


end % end of main function