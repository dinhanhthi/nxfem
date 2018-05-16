function RHS = getLvChopp07(tris,CT,msh,pa,uold,wS)
% load vector for v: - int_Omg bet*u*phi
% in this case w=-bet*uold, g(u)=1
% related file: main_chopp2007.m
% State: 
% Input: 
% Output: global load vector F (nNodes+nNew) x 1

CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
newNodes = msh.newNodes;
typeG = 5; % g(u)=1, cf. defG.m
% REMEMBER to change typeG in getGMgPP (or other global matrix) also

%% =======================================================================
% GET COMPONENTS
% ========================================================================

%-------------------------------------------------------------------------
% There is no term f(x,y)*phi (f=0 in this case)
%-------------------------------------------------------------------------

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
i1 = ig1; i2 = ig2; ic = igc;
f1 = fg1; f2 = fg2; fc1 = fgc1; fc2 = fgc2;



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