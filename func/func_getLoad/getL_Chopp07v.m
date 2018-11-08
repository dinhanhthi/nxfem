function RHS = getL_Chopp07v(msh,pa,tris,CT,uold)
% Load vector for v in chopp07
% SHAPE: - int_Omg bet*u*phi
% Based on getLf, change in getPf
% Related file: main_chopp2007.m, getGM_chopp07v.m
% State: (old) checked with old file before rewrite
% Input: - uold : solution uh, 4 component got from getWsep, including -beta
% Output: global load vector F (nNodes+nNew) x 1

CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
newNodes = msh.newNodes;

%% =======================================================================
% GET COMPONENTS
% ========================================================================

%-------------------------------------------------------------------------
% Term int -bet*u*phi
%-------------------------------------------------------------------------
% get P
sol.u = uold; % from getWsep, uold here already includes "-beta*"
func.gu = @(u,pa) u; % g(u) = u
P = getPf(msh,pa,tris,CT,sol,func);
% get indices
[i1,f1] = getfPhiNCTs(msh,pa,NCTs1,P.NC1); % NCTs1
[i2,f2] = getfPhiNCTs(msh,pa,NCTs2,P.NC2); % NCTs2
[ic,fc1,fc2] = getfPhiCTs(msh,pa,CTs,CT,P); % CTs



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
