function mLG2 = getMatrixL2Grad(omg1NCTs,omg2NCTs,nodesOmg2GamCTs,...
                    CTs,areaChildCTs,areaCTs,msh,pa)
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
[i1,j1,v1] = getMatrixGradGradNCTs(omg1NCTs,k1,msh);
[i2,j2,v2] = getMatrixGradGradNCTs(omg2NCTs,k2,msh);
[it,jt,vt1,vt2] = getMatrixGradGradCTs(CTs,areaChildCTs,areaCTs,k1,k2,msh);



%% ========================================================
% BUID MATRIX
% =========================================================
ii = [i1;it]; jj = [j1;jt]; vv = [v1;vt1]; % omg1NCTs & CTs1
itmp = newNodes(it); jtmp = newNodes(jt);
ii = [ii;itmp]; jj = [jj;jtmp]; vv = [vv;vt2]; % CTs2
% Replace all nodes in i2 and j2 which are also in nodesOmg2GamCTs
%   with the new values
tmp = ismember(i2,nodesOmg2GamCTs); % column-array
i2(tmp) = newNodes(i2(tmp)); % column-array
tmp = ismember(j2,nodesOmg2GamCTs); % column-array
j2(tmp) = newNodes(j2(tmp)); % column-array
ii = [ii;i2]; jj = [jj;j2]; vv = [vv;v2]; % omg2NCTs

mLG2 = sparse(ii,jj,vv);

end