function globalMatrix = getGlobalMatrix(omg1NCTs,omg2NCTs,CTs,phi,hT,...
       areaChildCTs,areaCTs,iPs,uNormalsCT,nodesCTsInOmg2OnGam,msh,pa,cp)
% Find the global stiffness matrix
% Input: - triangles on not cut triangles
%        - information about cut triangles
%	 	 - nodesOmg2GamCTs : to get new nodes in CTs2
% Output: global stiffness matrix A

newNodes = msh.newNodes; % convert i to k(i)

%% ========================================================
% GET MATRICES
% =========================================================
[mNCTs1,mNCTs2] = getMatrixNCTs(omg1NCTs,omg2NCTs,msh,pa); % on NCTs
[mCTs1,mCTs2,mCTs3,mCTs4] = getMatrixCTs(CTs,areaChildCTs,...
                            areaCTs,iPs,uNormalsCT,msh,pa,cp); % on CTs


%% ========================================================
% GET TRIPLETS
% They are column-arrays
% =========================================================
[i1,j1,v1] = find(mNCTs1); [i2,j2,v2] = find(mNCTs2); % on NCTs
% [i1,j1,v1] = getMatrixGradGradNCTs(omg1NCTs,pa.kk1,msh); % NCTs1
% [i2,j2,v2] = getMatrixGradGradNCTs(omg2NCTs,pa.kk2,msh); % NCTs2

[it1,jt1,vt1] = find(mCTs1); [it2,jt2,vt2] = find(mCTs2); % on CTs
[it3,jt3,vt3] = find(mCTs3); [it4,jt4,vt4] = find(mCTs4); % on CTs



%% ========================================================
% Get from basis whose support are on Omg1 (mNCTs1,mCTs1)
% We keep all indices
% =========================================================
iG = [i1;it1]; jG = [j1;jt1]; vG = [v1;vt1];



%% ========================================================
% Get from basis whose support are on Omg2 (mNCTs2,mCTs2)
% =========================================================

% Replace all nodes in i2 and j2 which are also in nodesCTsInOmg2OnGam with
%   the new values
tmp = ismember(i2,nodesCTsInOmg2OnGam); % column-array
i2(tmp) = newNodes(i2(tmp)); % column-array
tmp = ismember(j2,nodesCTsInOmg2OnGam); % column-array
j2(tmp) = newNodes(j2(tmp)); % column-array
% Add to global matrix
iG = [iG;i2]; jG = [jG;j2]; vG = [vG;v2];

% Replace all nodes in it2, jt2 with the new values 
% (because these nodes are in the cut triangles region)
it2 = newNodes(it2); % column-array
jt2 = newNodes(jt2); % column-array
% Add to global matrix
iG = [iG;it2]; jG = [jG;jt2]; vG = [vG;vt2];



%% ========================================================
% Get from nodes which locate only on cut triangles (mCTs3,mCTs4)
% =========================================================

% For mCTs3, this is term A_{k(i)j}, we replace all nodes in it with the 
%   new values, keep jt as original ones
it3 = newNodes(it3); % column-array
% Add to global matrix
iG = [iG;it3]; jG = [jG;jt3]; vG = [vG;vt3];

% For mCTs4, this is term A_{ik(j)}, we replace all nodes in jt with the 
%   new values, keep it as original ones
jt4 = newNodes(jt4); % column-array
% Add to global matrix
iG = [iG;it4]; jG = [jG;jt4]; vG = [vG;vt4];


%% ========================================================
% Ghost penalty terms
% =========================================================
if pa.useGP
    [iGP,jGP,vGP] = getGhostPenalty(CTs,phi,hT,msh,pa);
    iG = [iG;iGP]; jG = [jG;jGP]; vG = [vG;vGP];
end


%% ========================================================
% GLOBAL MATRIX
% =========================================================
globalMatrix = sparse(iG,jG,vG);

end