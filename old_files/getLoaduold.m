function RHS = getLoaduold(NCTs1,NCTs2,CTs,iPs,typeCTs,nodeCTs,msh,pa,defF,uold,wS)
% Get the right hand side load vector for equation u
% State:
% Input: - load vector on not-cut triangles in each subdomain : FNC1, FNC2
%        - load vector on cut-triangles for 2 cases: FCT1,FCT2
%        - idx of nodes (around the interface) in each subdomain and on the interface
% Output: global load vector F (nNodes+nNew)x1

newNodes = msh.newNodes;

%% ========================================================
% GET COMPONENTS
% =========================================================
FNC1 = getLoaduNCTs(NCTs1,1,msh,pa,defF,uold.omg1,wS.omg1); % NCTs1
FNC2 = getLoaduNCTs(NCTs2,2,msh,pa,defF,uold.omg2,wS.omg2); % NCTs2
[FCT1,FCT2] = getLoaduCTs(CTs,iPs,typeCTs,nodeCTs,msh,pa,defF,uold,wS); % CTs

%% ========================================================
% NOT CUT TRIANGLES in OMG1
% =========================================================
nFNC1 = size(FNC1,1);
ii = 1:nFNC1; % row array
ff = FNC1'; % because FNC1, FNC2 are column-arrays



%% ========================================================
% NOT CUT TRIANGLES in OMG2
% =========================================================
nFNC2 = size(FNC2,1);
itmp = 1:nFNC2;
% Get nodes in Omg2 and on Interface of the cut triangles
nodeOmg2IntCTs = [nodeCTs.Omg2,nodeCTs.onG]; % row-array
% Replace all nodes in itmp which are also in nodeOmg2Int with the new values
tmp = ismember(itmp,nodeOmg2IntCTs);
itmp(tmp) = newNodes(itmp(tmp)); % row-array
ii = [ii,itmp];
ff = [ff,FNC2']; % because FNC1, FNC2 are column-arrays



%% ========================================================
% CUT TRIANGLES
% =========================================================
% in Omg1
nFCT1 = size(FCT1,1);
ii = [ii,1:nFCT1];
ff = [ff,FCT1'];
% in Omg2
nFCT2 = size(FCT2,1);
itmp2 = 1:nFCT2;
tmp = ismember(itmp2,nodeCTs.all); % find nodes around interface
itmp2(tmp) = newNodes(itmp2(tmp)); % row-array
ii = [ii,itmp2];
ff = [ff,FCT2'];



%% ========================================================
% GLOBAL LOAD VECTOR
% =========================================================
ii = ii'; % transform to column array
RHS = accumarray(ii,ff); % column array

end