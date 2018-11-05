function RHS = getLfgu(msh,pa,tris,CT,uold,defF,defGu,coef1,coef2)
% Load vector for v in chopp07
% SHAPE: int_Omg f(u)*phi + int_Omg coef*g(u)*phi
% Based on getL_Chopp07v and getLf
% Related file: main_chopp06combine.m (first used)
% State:
% Input: - uold : solution uh, 4 component got from getWsep, including -beta
%        - defGu(u,pa): function handle from defGu.m
%        - coef1, coef2 w.r.t Omg1 and Omg2 (including the sign)
%        - RHS f(x,y) from model_...
% Output: global load vector F (nNodes+nNew) x 1

CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
newNodes = msh.newNodes;

%% =======================================================================
% GET COMPONENTS
% ========================================================================

%-------------------------------------------------------------------------
% Term int f*phi
%-------------------------------------------------------------------------
% get P
sol = []; func.h = defF;
P = getPf(msh,pa,tris,CT,sol,func);
% get indices
[i1,f1] = getfPhiNCTs(msh,pa,NCTs1,P.NC1); % NCTs1
[i2,f2] = getfPhiNCTs(msh,pa,NCTs2,P.NC2); % NCTs2
[ic,fc1,fc2] = getfPhiCTs(msh,pa,CTs,CT,P); % CTs


%-------------------------------------------------------------------------
% Term int coef*g(u)*phi (sign goes with coef)
%-------------------------------------------------------------------------
% get P
sol.u = uold; % from getWsep, uold here already includes coef
func.gu = defGu;
func.h = @(x,y,pa,sub) findDefH(x,y,pa,sub,coef1,coef2);
P = getPf(msh,pa,tris,CT,sol,func);
% get indices
[ig1,fg1] = getfPhiNCTs(msh,pa,NCTs1,P.NC1); % NCTs1
[ig2,fg2] = getfPhiNCTs(msh,pa,NCTs2,P.NC2); % NCTs2
[igc,fgc1,fgc2] = getfPhiCTs(msh,pa,CTs,CT,P); % CTs
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


    function val = findDefH(xx,yy,pa,sub,c1,c2)
       if sub==1
           val = c1;
       else
           val = c2;
       end
    end

end
