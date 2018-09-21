function RHS = getLwgP(NCTs1,NCTs2,CTs,iPs,typeCTs,nodeCTs,msh,pa,...
                            defF,uold,wS)
% load vector for u (there's more -w*g(u)*phi)
% State: checked, alike "old" getLoadu
% Input: - load vector on not-cut triangles in each subdomain : FNC1, FNC2
%        - load vector on cut-triangles for 2 cases: FCT1,FCT2
%        - idx of nodes (around the interface) in each subdomain and on the interface
% Output: global load vector F (nNodes+nNew)x1

newNodes = msh.newNodes;

%% =======================================================================
% GET COMPONENTS
% ========================================================================

%-------------------------------------------------------------------------
% Term int f(x,y)*phi
%-------------------------------------------------------------------------
[i1,f1] = getLFphiNCTs(NCTs1,msh,pa,defF,1); % NCTs1
[i2,f2] = getLFphiNCTs(NCTs2,msh,pa,defF,2); % NCTs2
[ic,fc1,fc2] = getLFphiCTs(CTs,iPs,typeCTs,nodeCTs,msh,pa,defF); % CTs

%-------------------------------------------------------------------------
% Term int wg(u)*phi
%-------------------------------------------------------------------------
[ig1,fg1] = getLwgNCTs(NCTs1,msh,pa,uold.omg1,wS.omg1,1); % NCTs1
[ig2,fg2] = getLwgNCTs(NCTs2,msh,pa,uold.omg2,wS.omg2,1); % NCTs2
[igc,fgc1,fgc2] = getLwgCTs(CTs,iPs,typeCTs,nodeCTs,msh,pa,uold,wS,1); % CTs
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