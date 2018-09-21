function [ii,jj,vv] = getTriplegPPNCTs(NCTs,uold,kk,msh,pa,typeG)
% Get triplets of int_NCTs(kk*g(uold)*phi*phi)
% Input: - information of NCTs of a specific subdomain
%        - coefficients kk on this subdomain
%        - uold
%        - type of g(u), see file defG.m
% Output: triplet on NCTs wrt each subdomain (column-arrays)

nNCTs = size(NCTs,2); % number of all not cut triangles in Omg1
ii = zeros(9*nNCTs,1); jj = zeros(9*nNCTs,1); vv = zeros(9*nNCTs,1);
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D

idx=1;
for t=1:nNCTs
    triangle = NCTs(:,t);
    uoldT(1:3) = uold(NCTs(1:3,t)); % uold's values at vertices of triangle
    for i=1:3
        for j=1:3
            ii(idx) = NCTs(i,t);
            jj(idx) = NCTs(j,t);
            vv(idx) = getTriplegPPWhole(triangle,uoldT,j,i,...
                            dim,deg,msh,typeG,kk);
            idx = idx+1;
        end
    end
end

end