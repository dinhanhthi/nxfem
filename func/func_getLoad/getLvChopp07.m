function RHS = getLvChopp07(msh,pa,tris,CT,wold)
% NEED TO BE EDITED W.R.T THE NEW REWRITE!!!
% Load vector for v: - int_Omg bet*u*phi
% We apply the idea from getLfwg but in this case w=-bet*uold, g(u)=1
%   (because we cannot defind g(u) with beta)
% Related file: main_chopp2007.m
% State: checked with old file before rewrite
% Input: 
% Output: global load vector F (nNodes+nNew) x 1

CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
newNodes = msh.newNodes;

%% =======================================================================
% GET COMPONENTS
% ========================================================================

%-------------------------------------------------------------------------
% There is no term f(x,y)*phi (f=0 in this case)
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% Term int wg(u)*phi
% w=bet*uold, g(u)=1
%-------------------------------------------------------------------------
defFw = @(w) w; % f(w) = w
[ig1,fg1] = getfPhiNCTs(msh,pa,NCTs1,'w',wold.omg1,'fw',defFw); % NCTs1
[ig2,fg2] = getfPhiNCTs(msh,pa,NCTs2,'w',wold.omg2,'fw',defFw); % NCTs2
[igc1,fgc1] = getfPhiCTs(msh,pa,CTs,CT,1,'w',wold.ct1,'fw',defFw); % CTs1
[igc2,fgc2] = getfPhiCTs(msh,pa,CTs,CT,2,'w',wold.ct2,'fw',defFw); % CTs2

% sign "-" in load vector form
fg1=-fg1; fg2=-fg2; fgc1=-fgc1; fgc2=-fgc2;

% add to terms f*phi
% (in the future, if there are more terms, just do like this)
i1 = ig1; i2 = ig2; ic1 = igc1; ic2 = igc2;
f1 = fg1; f2 = fg2; fc1 = fgc1; fc2 = fgc2;



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
