function segment = getMeshSegmentOnGh(msh,pa,phi)
% Get all segment of mesh whose endpoints located on Gamh
% Input: 
% Output: vector contains all segment (2 idx point x number of segment)

phi(abs(phi)<pa.tol)=0; % find phi which are very small (~0) and set to 0

phiTris = phi(msh.t(1:3,:)); % phi's values in size of msh.t
eq0idx3 = ceil(find(abs(phiTris)<pa.tol)/3);
[~,tmp] = unique(eq0idx3);
tmp = setdiff(1:numel(eq0idx3),tmp);
dupPointInTri = eq0idx3(tmp); % idx of triangle containing segment of Gamh (column array)
segment = zeros(2,numel(dupPointInTri));
for i=1:numel(dupPointInTri)
    segment(:,i) = msh.t(abs(phiTris(:,dupPointInTri(i)))<pa.tol,dupPointInTri(i));
end
segment = sort(segment);
segment = unique(segment','rows');
segment = segment';


end