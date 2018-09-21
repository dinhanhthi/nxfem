function mLG2 = getMatrixL2Grad(NCTs1,NCTs2,CTs,areaChildCTs,msh,pa)
% Find semi-norm |k^{1/2}grad|_L2
% State: checked
% Input: - information of NCTs and CTs
% Output: matrix to be used in err = E*A*E'


%% ========================================================
% INFORMATION
% =========================================================
k1=pa.kk1; k2=pa.kk2;
% k1=1; k2=1;
newNodes = msh.newNodes;



%% ========================================================
% GET TRIPLETS
% all are column-arrays
% =========================================================
[i1,j1,v1] = getTripleGGNCTs(NCTs1,k1,msh);
[i2,j2,v2] = getTripleGGNCTs(NCTs2,k2,msh);
[it,jt,vt1,vt2] = getTripleGGCTs(CTs,areaChildCTs,k1,k2,msh);



%% ========================================================
% BUID MATRIX
% =========================================================
ii = [i1;it]; jj = [j1;jt]; vv = [v1;vt1]; % omg1NCTs & CTs1
itmp = newNodes(it); jtmp = newNodes(jt);
ii = [ii;itmp]; jj = [jj;jtmp]; vv = [vv;vt2]; % CTs2
% Replace all nodes in i2 and j2 which are also in nodesOmg2GamCTs
%   with the new values
tmp = ismember(i2,msh.node.CT.omg2); % column-array
i2(tmp) = newNodes(i2(tmp)); % column-array
tmp = ismember(j2,msh.node.CT.omg2); % column-array
j2(tmp) = newNodes(j2(tmp)); % column-array
ii = [ii;i2]; jj = [jj;j2]; vv = [vv;v2]; % omg2NCTs

mLG2 = sparse(ii,jj,vv);

end