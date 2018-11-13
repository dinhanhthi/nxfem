function [ii,jj,vv] = getMatrixGijNCTs(NCTs,areaNCTs,typeG,pa)
% get triplet matrix of term int_Omg(G(phi_i,phi_j)) on not cut triangles
% Input: - information of NCTs of a specific subdomain
%        - type of function G (typeG)
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
    for i=1:3 % 3 vertices
       for j=1:3 % 3 vertices
           ii(idx) = NCTs(i,t);
           jj(idx) = NCTs(j,t);
           vv(idx) = getGij_WholeTri(areaNCTs(1,t),j,i,typeG,dim,deg); % A_ij=A(i,j)
           idx = idx+1;
       end % end for i
    end % end for j
end % end for nNCTs

end