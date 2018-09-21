function metaK = getMatEstEtaK(NCTs1,NCTs2,iPs,CTs,typeCTs,nodeCTs,pa,msh)
% Try to find a suitable value for lambda by using the of posteriori error
%   estimator. The form of estimator is given in page 33, Barrau's thesis.
% IMPORTANT: NEED TO CHANGE TO MORE GENERAL CASE (don't use Gij anymore),
% do this later when we need etak again
% This function computes the matrix to find the value of term etaK
% The idea is the same with function getMatrixL2 but we cannot use this
%   function in this case
% Input: 
% Output: matrix to be used in finding value of etaK

newNodes = msh.newNodes; % convert i to k(i)
hT = msh.hT;

%% ========================================================
% GET TRIPLETS
% all are column-arrays
% =========================================================
[i1,j1,v1] = getMatEstEtaKNCTs(NCTs1,hT,pa.kk1,pa,msh); % in Omg1
[i2,j2,v2] = getMatEstEtaKNCTs(NCTs2,hT,pa.kk2,pa,msh); % in Omg2
[it,jt,vt1,vt2] = getMatEstEtaKCTs(CTs,iPs,typeCTs,nodeCTs,hT,pa,msh);% CTs


                                    
%% ========================================================
% BUID MATRIX
% =========================================================
ii = [i1;it]; jj = [j1;jt]; vv = [v1;vt1]; % omg1NCTs & CTs1
itmp = newNodes(it); jtmp = newNodes(jt);
ii = [ii;itmp]; jj = [jj;jtmp]; vv = [vv;vt2]; % CTs2
% Replace all nodes in i2 and j2 which are also in nodesCTsInOmg2OnGam
%      with the new values
tmp = ismember(i2,msh.node.CT.omg2); % column-array
i2(tmp) = newNodes(i2(tmp)); % column-array
tmp = ismember(j2,msh.node.CT.omg2); % column-array
j2(tmp) = newNodes(j2(tmp)); % column-array
ii = [ii;i2]; jj = [jj;j2]; vv = [vv;v2]; % omg2NCTs

metaK = sparse(ii,jj,vv);
end