function RHS = getLfwg(tris,CT,msh,pa,defF,uold,wS)
% load vector for u: int_Omg f*phi - int_Omg wg(u)*phi
% related file: getGMgPP, main_sys_linda
% State: checked with getLwgP (old file), alike upto e-14
% Input:
% Output: global load vector F (nNodes+nNew)x1

CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
newNodes = msh.newNodes;
typeG = 1; % g(u) as normal, cf. defG.m
% REMEMBER to change typeG in getGMgPP (or other global matrix) also

%% =======================================================================
% GET COMPONENTS
% ========================================================================

%-------------------------------------------------------------------------
% Term int f(x,y)*phi
%-------------------------------------------------------------------------
Pf = getPf(tris,CT,msh,pa,defF);
[i1,f1] = getLNCTs(NCTs1,msh,pa,Pf.NC1); % NCTs1
[i2,f2] = getLNCTs(NCTs2,msh,pa,Pf.NC2); % NCTs2
[ic,fc1,fc2] = getLCTs(CTs,CT,msh,pa,Pf); % CTs

%-------------------------------------------------------------------------
% Term int wg(u)*phi
%-------------------------------------------------------------------------
Pwg = getPwg(tris,CT,uold,wS,msh,pa,typeG);
[ig1,fg1] = getLNCTs(NCTs1,msh,pa,Pwg.NC1); % NCTs1
[ig2,fg2] = getLNCTs(NCTs2,msh,pa,Pwg.NC2); % NCTs2
[igc,fgc1,fgc2] = getLCTs(CTs,CT,msh,pa,Pwg); % CTs
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