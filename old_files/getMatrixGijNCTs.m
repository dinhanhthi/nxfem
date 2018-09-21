function [ii,jj,vv] = getMatrixGijNCTs(NCTs,pa,msh)
% get triplet matrix of term int_Omg(G(phi_i,phi_j)) on not cut triangles
% Input: - information of NCTs of a specific subdomain
% Output: triplet in NCTs wrt this subdomain (column-arrays)


%% ========================================================
% INFORMATION
% =========================================================
nNCTs = size(NCTs,2); % number of all not cut triangles
ii = zeros(9*nNCTs,1); % column-array
jj = zeros(9*nNCTs,1); % column-array
vv = zeros(9*nNCTs,1); % column-array


%% ========================================================
% GET TRIPLET
% =========================================================
idx = 1;
for t=1:nNCTs
    triangle = NCTs(:,t);
    for i=1:3 % 3 vertices
       for j=1:3 % 3 vertices
           ii(idx) = NCTs(i,t);
           jj(idx) = NCTs(j,t);
           vv(idx) = getGij_WholeTri(triangle,j,i,pa,msh); % A_ij=A(i,j)
           idx = idx+1;
       end % end for i
    end % end for j
end % end for nNCTs

end