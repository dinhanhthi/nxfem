function RHS = getLfwg(msh,pa,tris,CT,uold,wold,defF,defG)
% load vector for u: int_Omg f*phi - int_Omg wg(u)*phi
% Related file: getGMgPP, main_sys_linda
% Old file: getLwg_old.m
% State: - checked wrt old file getLfwg.m
%        - checked with getLwgP (old file), alike upto e-14
% Input: - mesh msh and all triangles tris
%        - CTs' components
%        - uold, wold has .omg1, .omg2, .ct1, .ct2 (cf. getWsep)
%        - defF(x,y,pa,sub)
%        - defG: g(u)
% Output: global load vector F (nNodes+nNew)x1

CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
newNodes = msh.newNodes;


%% =======================================================================
% GET COMPONENTS
% ========================================================================

defH1 = @(x,y,pa) defF(x,y,pa,1);
defH2 = @(x,y,pa) defF(x,y,pa,2);

%-------------------------------------------------------------------------
% Term int f(x,y)*phi
%-------------------------------------------------------------------------
[i1,f1] = getfPhiNCTs(msh,pa,NCTs1,'h',defH1); % NCTs1
[i2,f2] = getfPhiNCTs(msh,pa,NCTs2,'h',defH2); % NCTs2
[ic1,fc1] = getfPhiCTs(msh,pa,CTs,CT,1,'h',defH1); % CTs1
[ic2,fc2] = getfPhiCTs(msh,pa,CTs,CT,2,'h',defH2); % CTs2

%-------------------------------------------------------------------------
% Term int wg(u)*phi
%-------------------------------------------------------------------------
defFw = @(w) w; % f(w) = w
defGu = defG.change;
[ig1,fg1] = getfPhiNCTs(msh,pa,NCTs1,'u',uold.omg1,'gu',defGu,...
                'w',wold.omg1,'fw',defFw); % NCTs1
[ig2,fg2] = getfPhiNCTs(msh,pa,NCTs2,'u',uold.omg2,'gu',defGu,...
                'w',wold.omg2,'fw',defFw); % NCTs2
[igc1,fgc1] = getfPhiCTs(msh,pa,CTs,CT,1,'u',uold.ct1,'gu',defGu,...
                'w',wold.ct1,'fw',defFw); % CTs1
[igc2,fgc2] = getfPhiCTs(msh,pa,CTs,CT,2,'u',uold.ct2,'gu',defGu,...
                'w',wold.ct2,'fw',defFw); % CTs2

% sign "-" in load vector form
fg1=-fg1; fg2=-fg2; fgc1=-fgc1; fgc2=-fgc2;

% add to terms f*phi
% (in the future, if there are more terms, just do like this)
i1 = [i1;ig1]; i2 = [i2;ig2]; ic1 = [ic1;igc1]; ic2 = [ic2;igc2];
f1 = [f1;fg1]; f2 = [f2;fg2]; fc1 = [fc1;fgc1]; fc2 = [fc2;fgc2];



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