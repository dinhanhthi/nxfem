function [ii,ff] = getLwgNCTs(NCTs,msh,pa,uold,wS,typeG)
% Get int_NCTs wg(u)*phi for equation of u
% State: 
% Input: - NCTs
%        - wS, uold (in seperated subdomain)
% Output: ii (nodes), ff (values at nodes). column arrays

nNCTs = size(NCTs,2); % number of triangleNC

ii = zeros(3*nNCTs,1); ff = zeros(3*nNCTs,1);
idx=1;

for t=1:nNCTs
    uoldT(1:3) = uold(NCTs(1:3,t)); % uold's values at vertices of triangle
    wST(1:3) = wS(NCTs(1:3,t));
    tri = NCTs(:,t);
    for i=1:3 % 3 vertices
        ii(idx) = NCTs(i,t);
        ff(idx) = getLwgWhole(tri,i,pa,msh,uoldT,wST,typeG);
        idx = idx+1;
    end % end for vertices
end % end for nTriNC

end