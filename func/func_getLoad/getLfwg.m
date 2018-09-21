function RHS = getLfwg(msh,pa,tris,CT,uold,wold,defF,defG)
% Load vector for u: int_Omg f(x,y)*phi - int_Omg wg(u)*phi
% Note: rewriten after 1st rewrite with inputParser (worst performance)
% related file: getGMgPP, main_sys_linda
% State: - checked, faster then the old inputParser
%        - (old) checked with getLwgP (old file), alike upto e-14
% Input: - tris : all triangles
%        - CT: all components of cut triangles
%        - uold, wold: .omg1, .omg2, .ct1, .ct2, all are in stdFEM
%        - defG: g(u)
%        - defF(x,y,pa,sub) : function handle
% Output: global load vector F (nNodes+nNew)x1

CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
newNodes = msh.newNodes;


%% =======================================================================
% GET COMPONENTS
% ========================================================================

%-------------------------------------------------------------------------
% Term int f(x,y)*phi
%-------------------------------------------------------------------------
% get P like in getLf.m
sol = []; func.h = defF;
Pf = getPf(msh,pa,tris,CT,sol,func);
[i1,f1] = getfPhiNCTs(msh,pa,NCTs1,Pf.NC1); % NCTs1
[i2,f2] = getfPhiNCTs(msh,pa,NCTs2,Pf.NC2); % NCTs2
[ic,fc1,fc2] = getfPhiCTs(msh,pa,CTs,CT,Pf); % CTs

%-------------------------------------------------------------------------
% Term int wg(u)*phi
%-------------------------------------------------------------------------
sol.u = uold; sol.w = wold; 
func2.gu = defG.change;
func2.fw = @(w) w;
Pwg = getPf(msh,pa,tris,CT,sol,func2);
[ig1,fg1] = getfPhiNCTs(msh,pa,NCTs1,Pwg.NC1); % NCTs1
[ig2,fg2] = getfPhiNCTs(msh,pa,NCTs2,Pwg.NC2); % NCTs2
[igc,fgc1,fgc2] = getfPhiCTs(msh,pa,CTs,CT,Pwg); % CTs
% sign "-" in load vector form
fg1=-fg1; fg2=-fg2; fgc1=-fgc1; fgc2=-fgc2;
% add to terms f*phi
% (in the future, if there are more terms, just do like this)
i1 = [i1;ig1]; i2 = [i2;ig2]; ic = [ic;igc];
f1 = [f1;fg1]; f2 = [f2;fg2]; fc1 = [fc1;fgc1]; fc2 = [fc2;fgc2];



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
