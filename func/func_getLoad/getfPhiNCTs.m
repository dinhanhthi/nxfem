function [ii,ff] = getfPhiNCTs(msh,pa,NCTs,P)
% int_Omg P*phi on NCTs (load vector)
% This file is used for general P (the size of P depends on number of
%       Gaussian points (pa.degN)
% Note: rewriten after 1st rewrite with inputParser (worst performance)
% State: checked with getLf
% Related: getPfPhi
% Input: - not cut triangles NCTs
%        - component P: nT x nwt
% Output: ii (nodes), ff (values at nodes). column arrays

nNCTs = size(NCTs,2); % number of triangleNC

% default P
if isempty(P)
    dim=2; deg=pa.degN; 
    [wt,~] = getGaussQuad(dim,deg); 
    nwt = size(wt,2);
    P = ones(nNCTs,nwt);
end

ii = zeros(3*nNCTs,1); ff = zeros(3*nNCTs,1);
idx=1;
for t=1:nNCTs
    triangle = NCTs(:,t);
    PP = P(t,:);
    for i=1:3 % 3 vertices
        ii(idx) = NCTs(i,t);
        ff(idx) = getfPhiWhole(msh,pa,triangle,i,PP);
        idx = idx+1;
    end % end for vertices
end % end for nTriNC
end
