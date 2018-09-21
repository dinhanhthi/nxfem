
function RHS = getLf(tris,CT,msh,pa,defF)
% Get the load vector int_Omg f(x,y)*phi
% Related files: main_sys_linda, main_eachStep,...
% State: checked with getLfPhi (old file), alike upto e-16
% Input: - load vector on not-cut triangles in each subdomain : FNC1, FNC2
%        - load vector on cut-triangles for 2 cases: FCT1,FCT2
%        - idx of nodes (around the interface) in each subdomain and on the interface
% Output: global load vector F (nNodes+nNew)x1

CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
newNodes = msh.newNodes;

%% =======================================================================
% GET COMPONENTS
% ========================================================================

%-------------------------------------------------------------------------
% Term int f*phi
%-------------------------------------------------------------------------
% get P
P = getPf(tris,CT,msh,pa,defF);
% get indices
[i1,f1] = getLNCTs(NCTs1,msh,pa,P.NC1); % NCTs1
[i2,f2] = getLNCTs(NCTs2,msh,pa,P.NC2); % NCTs2
[ic,fc1,fc2] = getLCTs(CTs,CT,msh,pa,P); % CTs



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
