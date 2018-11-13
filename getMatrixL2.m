function mL2 = getMatrixL2(omg1NCTs,omg2NCTs,areaNCTs1,areaNCTs2,...
                iPs,CTs,typeCTs,nodeCTs,areaCTs,nodesOmg2GamCTs,msh,pa)
% Get matrix AL2 to find L2 norm error = int_Omg(phi_i*phi_j)
% State: checked
% Input: - all information abt not cut triangles
%        - all information abt cut triangles
%	 	 - nodesOmg2GamCTs : to get new nodes in CTs2
% Output: matrix to be used in errL2 = E*AL2*E'


%% ========================================================
% INFORMATION
% =========================================================
typeG = 1; % G(X,Y)=X*Y
newNodes=msh.newNodes;


%% ========================================================
% GET TRIPLETS
% all are column-arrays
% =========================================================
[i1,j1,v1] = getMatrixGijNCTs(omg1NCTs,areaNCTs1,typeG,pa); % in Omg1
[i2,j2,v2] = getMatrixGijNCTs(omg2NCTs,areaNCTs2,typeG,pa); % in Omg2
[it,jt,vt1,vt2] = getMatrixGijCTs(CTs,iPs,typeCTs,nodeCTs,...
                                        areaCTs,typeG,msh,pa); % in CTs

                                    

%% ========================================================
% BUID MATRIX
% =========================================================
ii = [i1;it]; jj = [j1;jt]; vv = [v1;vt1]; % omg1NCTs & CTs1
itmp = newNodes(it); jtmp = newNodes(jt);
ii = [ii;itmp]; jj = [jj;jtmp]; vv = [vv;vt2]; % CTs2
% Replace all nodes in i2 and j2 which are also in nodesOmg2GamCTs
%      with the new values
tmp = ismember(i2,nodesOmg2GamCTs); % column-array
i2(tmp) = newNodes(i2(tmp)); % column-array
tmp = ismember(j2,nodesOmg2GamCTs); % column-array
j2(tmp) = newNodes(j2(tmp)); % column-array
ii = [ii;i2]; jj = [jj;j2]; vv = [vv;v2]; % omg2NCTs

mL2 = sparse(ii,jj,vv);

end