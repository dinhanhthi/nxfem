function [ii,jj,vv] = getTriplePPNCTs(NCTs,msh,pa,K)
% Get triplets of int_NCTs(K*phi*phi)
% This file is used for general K (the size of K depends on number of
%       Gaussian points (pa.degN)
% Input: - information of NCTs of a specific subdomain
%        - component K
% Output: triplet on NCTs wrt each subdomain (column-arrays)

nNCTs = size(NCTs,2); % number of all not cut triangles in Omg1
ii = zeros(9*nNCTs,1); jj = zeros(9*nNCTs,1); vv = zeros(9*nNCTs,1);

idx=1;
for t=1:nNCTs
    triangle = NCTs(:,t);
    if isempty(K)
        KK = [];
    else
        KK = K(t,:);
    end
    for i=1:3
        for j=1:3
            ii(idx) = NCTs(i,t);
            jj(idx) = NCTs(j,t);
            vv(idx) = getTriplePPWhole(triangle,j,i,pa,msh,KK);
            idx = idx+1;
        end
    end
end

end