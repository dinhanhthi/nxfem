function RHS = getLfphi(NCTs1,NCTs2,CTs,iPs,typeCTs,nodeCTs,msh,pa,defF)
% Get the right hand side load vector for equation tw (there is only f*phi)
% State: alike with "old" getLoadtw
% Input: - load vector on not-cut triangles in each subdomain : FNC1, FNC2
%        - load vector on cut-triangles for 2 cases: FCT1,FCT2
%        - idx of nodes (around the interface) in each subdomain and on the interface
% Output: global load vector F (nNodes+nNew)x1

newNodes = msh.newNodes;

%% =======================================================================
% GET COMPONENTS
% ========================================================================

%-------------------------------------------------------------------------
% Term int f*phi
%-------------------------------------------------------------------------
[i1,f1] = getLFphiNCTs(NCTs1,msh,pa,defF,1); % NCTs1
[i2,f2] = getLFphiNCTs(NCTs2,msh,pa,defF,2); % NCTs2
[ic,fc1,fc2] = getLFphiCTs(CTs,iPs,typeCTs,nodeCTs,msh,pa,defF); % CTs



%% =======================================================================
% BUILD ii and ff
% ========================================================================

%-------------------------------------------------------------------------
% Keep all ii and ff for NCTs1 and CTs1
%-------------------------------------------------------------------------
ii = [i1;ic]; ff = [f1;fc1];

%-------------------------------------------------------------------------
% NCTs2
%-------------------------------------------------------------------------
tmp = ismember(i2,msh.node.CT.omg2); % column-array
i2(tmp) = newNodes(i2(tmp)); % column-array
ii = [ii;i2]; ff = [ff;f2];

%-------------------------------------------------------------------------
% CTs2
%-------------------------------------------------------------------------
ic = newNodes(ic);
ii = [ii;ic]; ff = [ff;fc2];



%% ========================================================
% GLOBAL LOAD VECTOR
% =========================================================
RHS = accumarray(ii,ff); % column array

end