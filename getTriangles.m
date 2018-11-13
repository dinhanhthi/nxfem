function [CTs,omg1NCTs,omg2NCTs] = getTriangles(phi,msh,pa)
% Find the elements cut by zero level set function
% Input: triangles and points of the mesh AND level set function phi
% Output: - triangles cut by interface and their idx (% 5 x nTs)
%         - triangles not cut by interface and in Omg1 (% 5 x nTs)
%         - triangles not cut by interface and in Omg2 (% 5 x nTs)

%% ========================================================
% Preliminary
% =========================================================
triangles = msh.t;
phiOnTri = phi(triangles(1:3,:)); % phi on vertices of each triangle
maxPhiOnTri = max(phiOnTri); % max phi at vertice on each triangle
minPhiOnTri = min(phiOnTri); % min phi at vertice on each triangle
maxMin = maxPhiOnTri.*minPhiOnTri;


%% ========================================================
% triangles cut by interface
% =========================================================
% the last element of each triangle contains idx of this triangle
% idxCTs = find(maxMin<0);
idxCTs = find((maxMin<0)&(abs(maxPhiOnTri)>pa.tol)&(abs(minPhiOnTri)>pa.tol));
CTs = [triangles(:,idxCTs);idxCTs]; % 5 x nTs


%% ========================================================
% triangles not cut by interface and in Omg1
% ========================================================= 
% the last element of each triangle contains idx of this triangle
% idxOmg1NCTS = find(maxPhiOnTri<=0);
idxOmg1NCTS = find((maxPhiOnTri<0)|(abs(maxPhiOnTri)<pa.tol));
omg1NCTs = [triangles(:,idxOmg1NCTS);idxOmg1NCTS]; % 5 x nTs
% areaNCTs1 = getAreaTris(omg1NCTs,msh);


%% ========================================================
% triangles not cut by interface and in Omg2
% ========================================================= 
% the last element of each triangle contains idx of this triangle
% idxOmg2NCTs = find(minPhiOnTri>=0);
idxOmg2NCTs = find((minPhiOnTri>0)|(abs(minPhiOnTri)<pa.tol));
omg2NCTs = [triangles(:,idxOmg2NCTs);idxOmg2NCTs]; % 5 x nTs
% areaNCTs2 = getAreaTris(omg2NCTs,msh);

end
