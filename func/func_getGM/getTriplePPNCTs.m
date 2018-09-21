function [ii,jj,vv] = getTriplePPNCTs(msh,pa,NCTs,K)
% Get triplets of int_NCTs(K*phi*phi)
% This file is used for general K (the size of K depends on number of
%       Gaussian points (pa.degN)
% Note: rewriten after 1st rewrite with inputParser (worst performance)
% State: - checked with getGMgPP
% Input: - not cut triangles NCTs
%        - component K: nT x nwt
% Output: triplet on NCTs wrt each subdomain (column-arrays)

nNCTs = size(NCTs,2); % number of all not cut triangles in Omg1
ii = zeros(9*nNCTs,1); jj = zeros(9*nNCTs,1); vv = zeros(9*nNCTs,1);

% default K
if isempty(K)
    dim=2; deg=pa.degN; 
    [wt,~] = getGaussQuad(dim,deg); 
    nwt = size(wt,2);
    K = ones(nNCTs,nwt);
end


idx=1;
for t=1:nNCTs
    triangle = NCTs(:,t);
    KK = K(t,:);
    for i=1:3
        for j=1:3
            ii(idx) = NCTs(i,t);
            jj(idx) = NCTs(j,t);
            vv(idx) = getTriplePPWhole(msh,pa,triangle,j,i,KK);
            idx = idx+1;
        end
    end
end

end
