function [ii,jj,vv] = getMatEstEtaKNCTs(NCTs,hT,kk,pa,msh)
% get triplet on not cut triangles only
% it will be used in function getMatEstEtaK
% this file is diff from file getMatrixGijNCTs becaue it contains hT
% Input: - information of NCTs of a specific subdomain
%        - type of function G (typeG)
%        - diffusion coefficient kk
%        - diameter of all triangles in NCTs (vector hT)
% Output: triplet in NCTs wrt this subdomain (column-arrays)

%% ========================================================
% INFORMATION
% =========================================================
dim=2; deg=pa.degP2D; %Gaussian quadrature points in 2D
nNCTs = size(NCTs,2); % number of all not cut triangles
ii = zeros(9*nNCTs,1); % column-array
jj = zeros(9*nNCTs,1); % column-array
vv = zeros(9*nNCTs,1); % column-array

%% ========================================================
% GET TRIPLET
% =========================================================
idx = 1;
for t=1:nNCTs
    dK = hT(NCTs(5,t)); % diam of triangle K
    triangle = NCTs(:,t);
    for i=1:3 % 3 vertices
       for j=1:3 % 3 vertices
           ii(idx) = NCTs(i,t);
           jj(idx) = NCTs(j,t);
           % A_ij=A(i,j)
           vv(idx) = dK^2/kk*...
               getGij_WholeTri(triangle,j,i,dim,deg,msh);
           idx = idx+1;
       end % end for i
    end % end for j
end % end for nNCTs

end