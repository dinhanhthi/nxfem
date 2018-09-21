function FNC = getLoadNCTs(NCTs,areaNCTs,sub,msh,pa,defF)
% Get small-local load vector for triangles in Omega_sub and are 
%   not cut by interface
% NOTE that, this technique is GENERAL (doesn't depend on the number of
%               Gaussian points)
% State: checked by hand: code follows idea + check with double integral
% Input: - points and triangles in Omega_sub and are not cut by interface
%        - subdomain sub
%        - area of triangles in triangleNC
% Output: load vector defined on not-cut triangles in Omg_sub

% pa.degN: Gaussian quadrature points (for complicated function)
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
nNCTs = size(NCTs,2); % number of triangleNC

ii = zeros(3*nNCTs,1); ff = zeros(3*nNCTs,1);
idx=1;

for t=1:nNCTs
    for i=1:3 % 3 vertices
        ii(idx) = NCTs(i,t);
        ff(idx) = getLoadWholeTri(NCTs(:,t),...
                        sub,areaNCTs(1,t),i,dim,deg,msh,pa,defF);
        idx = idx+1;
    end % end for vertices
end % end for nTriNC

FNC = accumarray(ii,ff); % column array
end