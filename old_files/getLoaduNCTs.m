function FNC = getLoaduNCTs(NCTs,sub,msh,pa,defF,uold,wS)
% Get int_NCTs (fu-wg(u))phi for equation of u
% NOTE that, this technique is GENERAL (doesn't depend on the number of
%               Gaussian points)
% State:
% Input: - points and triangles in Omega_sub and are not cut by interface
%        - subdomain sub
%        - wS, uold (in seperated subdomain)
% Output: load vector defined on not-cut triangles in Omg_sub

% pa.degN: Gaussian quadrature points (for complicated function)
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
nNCTs = size(NCTs,2); % number of triangleNC

ii = zeros(3*nNCTs,1); ff = zeros(3*nNCTs,1);
idx=1;

for t=1:nNCTs
    uoldT(1:3) = uold(NCTs(1:3,t)); % uold's values at vertices of triangle
    wST(1:3) = wS(NCTs(1:3,t));
    tri = NCTs(:,t);
    for i=1:3 % 3 vertices
        ii(idx) = NCTs(i,t);
        ff(idx) = getLoaduWhole(tri,sub,i,dim,deg,msh,pa,defF,uoldT,wST);
        idx = idx+1;
    end % end for vertices
end % end for nTriNC

FNC = accumarray(ii,ff); % column array
end