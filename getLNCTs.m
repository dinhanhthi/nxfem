function [ii,ff] = getLNCTs(NCTs,msh,pa,P)
% int_Omg P*phi on NCTs (load vector)
% This file is used for general P (the size of P depends on number of
%       Gaussian points (pa.degN)
% State: alike getLoadCTs, diff at 16th digit after ","
% Input: - not cut triangles NCTs
%        - component P: nT x nwt
% Output: ii (nodes), ff (values at nodes). column arrays

nNCTs = size(NCTs,2); % number of triangleNC

ii = zeros(3*nNCTs,1); ff = zeros(3*nNCTs,1);
idx=1;
for t=1:nNCTs
    tri = NCTs(:,t);
    PP = P(t,:);
    for i=1:3 % 3 vertices
        ii(idx) = NCTs(i,t);
        ff(idx) = getLWhole(tri,i,msh,pa,PP);
        idx = idx+1;
    end % end for vertices
end % end for nTriNC
end