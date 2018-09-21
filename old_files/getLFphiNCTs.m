function [ii,ff] = getLFphiNCTs(NCTs,msh,pa,defF,kk)
% int_Omg f*phi on NCTs (load vector)
% State: alike getLoadCTs, diff at 16th digit after ","
% Input: - not cut triangles NCTs
%        - which domain? kk
%        - definition of f defF
% Output: ii (nodes), ff (values at nodes). column arrays

nNCTs = size(NCTs,2); % number of triangleNC

ii = zeros(3*nNCTs,1); ff = zeros(3*nNCTs,1);
idx=1;
for t=1:nNCTs
    triangle = NCTs(:,t);
    for i=1:3 % 3 vertices
        ii(idx) = NCTs(i,t);
        ff(idx) = getLFphiWhole(triangle,i,msh,pa,defF,kk);
        idx = idx+1;
    end % end for vertices
end % end for nTriNC
end