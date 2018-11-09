function val = getNormL2G(err,tris,areaChildCTs,msh,cp)
% Find semi-norm ||k^{1/2}grad||_L2
% if wanna take k=1, just modify cp before apply to getNormL2G
% State: checked
% Input: - information of NCTs and CTs
% Output: matrix to be used in err = E*A*E'

%-------------------------------------------------------------------------
% Get infor
%-------------------------------------------------------------------------
CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
newNodes = msh.newNodes;
if isempty(cp)
    k1=cp.kk1; k2=cp.kk2;
else
    k1=1; k2=1;
end

%-------------------------------------------------------------------------
% Get triplets
%-------------------------------------------------------------------------
[i1,j1,v1] = getTripleGGNCTs(NCTs1,k1,msh);
[i2,j2,v2] = getTripleGGNCTs(NCTs2,k2,msh);
[it,jt,vt1,vt2] = getTripleGGCTs(CTs,areaChildCTs,k1,k2,msh);

%-------------------------------------------------------------------------
% Build matrix AL2G
%-------------------------------------------------------------------------
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
mL2G = sparse(ii,jj,vv);

%-------------------------------------------------------------------------
% Get error L2G
%-------------------------------------------------------------------------
val = err'*mL2G*err;
val = sqrt(val);

end