function RHS = getLf(msh,pa,tris,CT,defF)
% Get the load vector int_Omg f(x,y)*phi
% Related files: main_sys_linda, main_eachStep,...
% State: - checked with getLf old
%        - checked with getLfPhi (old file), alike upto e-16
% Old file: getLf_old
% Input: - tris : all triangles
%        - CT: all components of cut triangles
%        - defF(x,y,pa,sub) : function handle
% Output: global load vector F (nNodes+nNew)x1

CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
newNodes = msh.newNodes;

%% =======================================================================
% GET COMPONENTS
% ========================================================================

defH1 = @(x,y,pa) defF(x,y,pa,1);
defH2 = @(x,y,pa) defF(x,y,pa,2);

%-------------------------------------------------------------------------
% Term int f*phi
%-------------------------------------------------------------------------
[i1,f1] = getfPhiNCTs(msh,pa,NCTs1,'h',defH1); % NCTs1
[i2,f2] = getfPhiNCTs(msh,pa,NCTs2,'h',defH2); % NCTs2
[ic1,fc1] = getfPhiCTs(msh,pa,CTs,CT,1,'h',defH1); % CTs1
[ic2,fc2] = getfPhiCTs(msh,pa,CTs,CT,2,'h',defH2); % CTs2



%% =======================================================================
% BUILD ii and ff
% ========================================================================

%-------------------------------------------------------------------------
% Keep all ii and ff for NCTs1 and CTs1
%-------------------------------------------------------------------------
ii = [i1;ic1]; ff = [f1;fc1];

%-------------------------------------------------------------------------
% NCTs2
%-------------------------------------------------------------------------
tmp = ismember(i2,msh.node.CT.omg2); % column-array
i2(tmp) = newNodes(i2(tmp)); % column-array
ii = [ii;i2]; ff = [ff;f2];

%-------------------------------------------------------------------------
% CTs2
%-------------------------------------------------------------------------
ic2 = newNodes(ic2);
ii = [ii;ic2]; ff = [ff;fc2];



%% ========================================================
% GLOBAL LOAD VECTOR
% =========================================================
RHS = accumarray(ii,ff); % column array


end
