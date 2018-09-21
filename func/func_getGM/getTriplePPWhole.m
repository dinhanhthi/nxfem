function val = getTriplePPWhole(msh,pa,triangle,i,j,K)
% Find int_NCTs(K*phi*phi) on whole triangle
% This file is used for general K (the size of K depends on number of
%       Gaussian points (dim,deg)
% Note: rewriten after 1st rewrite with inputParser (worst performance)
% state: - checked with getGMgPP
% Input: - which node i,j (1,2 or 3)
%        - which triangle? (with points to get vertices)
%        - component K: 1 x nwt
% Output: scalar value


%% Get info
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D (polynomial func)
[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points
points = msh.p;


%% default K
if isempty(K)
    K = ones(1,nwt);
end


%% find val
v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3
areaT = getAreaTri(v1,v2,v3); % area of triangle 

val=0;
for k=1:nwt
    [shFu,~,~] = getP1shapes(pt(1,k),pt(2,k)); % N_i at quadrature points
    val = val + K(k)*areaT*wt(k)*shFu(i)*shFu(j);
end


end
