function [ii,jj,vv] = getTripleGGNCTs(NCTs,kk,msh)
% Get triplets of int_NCTs(kk*gradPhi_i*gradPhi_j)
% Input: - information of NCTs of a specific subdomain
%        - coefficients kk on this subdomain
% Output: triplet on NCTs wrt each subdomain (column-arrays)


%% ========================================================
% PRELIMINARY
% =========================================================
nNCTs = size(NCTs,2); % number of all not cut triangles in Omg1
ii = zeros(9*nNCTs,1); jj = zeros(9*nNCTs,1); vv = zeros(9*nNCTs,1);
points=msh.p;

% normal vectors
% ve is a vector of 2 coordinates x nTriNC x 3 edges
ve(:,:,3) = points(:,NCTs(2,:)) - points(:,NCTs(1,:));
ve(:,:,1) = points(:,NCTs(3,:)) - points(:,NCTs(2,:));
ve(:,:,2) = points(:,NCTs(1,:)) - points(:,NCTs(3,:));
areaNCTs = 0.5*abs(ve(1,:,3).*ve(2,:,2) - ve(1,:,2).*ve(2,:,3));


%% ========================================================
% GET TRIPLET
% =========================================================
idx = 0;
for i=1:3
    for j=1:3
        % for Omg1NCTs
        ii(idx+1:idx+nNCTs) = NCTs(i,:);
        jj(idx+1:idx+nNCTs) = NCTs(j,:);
        % A_ij=A(j,i)
        vv(idx+1:idx+nNCTs) = kk*dot(ve(:,:,j),ve(:,:,i))./(4*areaNCTs); 
        idx = idx+nNCTs;
    end
end

end