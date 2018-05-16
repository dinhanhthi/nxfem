function vL2 = getNormL2(uh,tris,CT,msh,pa)
% Get the ||e||_L2 of a nxfem solution
% State: checked with the old method getMatrixL2.m
% Input: - solution in nxfem space: uh
%        - all needed info
% Output: value of L2 error

CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
newNodes = msh.newNodes;

%-------------------------------------------------------------------------
% Get K (= 1)
%-------------------------------------------------------------------------
nNCTs1 = size(NCTs1,2); nNCTs2 = size(NCTs2,2); nCTs = size(CTs,2);
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
[wt,~] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points
K.NC1 = ones(nNCTs1,nwt); K.NC2 = ones(nNCTs2,nwt);
K.CTW1 = ones(nCTs,nwt); K.CTW2 = K.CTW1; K.CTP1 = K.CTW1; K.CTP2 = K.CTW1;

%-------------------------------------------------------------------------
% Get triplets
%-------------------------------------------------------------------------
[i1,j1,v1] = getTriplePPNCTs(NCTs1,msh,pa,K.NC1); % NCTs1
[i2,j2,v2] = getTriplePPNCTs(NCTs2,msh,pa,K.NC2); % NCTs2
[it,jt,vt1,vt2] = getTriplePPCTs(CTs,CT,msh,pa,K); % CTs

%-------------------------------------------------------------------------
% Build matrix AL2
%-------------------------------------------------------------------------
ii = [i1;it]; jj = [j1;jt]; vv = [v1;vt1]; % NCTs2 & CTs1
itmp = newNodes(it); jtmp = newNodes(jt);
ii = [ii;itmp]; jj = [jj;jtmp]; vv = [vv;vt2]; % CTs2
% Replace all nodes in i2 and j2 which are also in msh.node.CT.omg2
%      with the new values
tmp = ismember(i2,msh.node.CT.omg2); % column-array
i2(tmp) = newNodes(i2(tmp)); % column-array
tmp = ismember(j2,msh.node.CT.omg2); % column-array
j2(tmp) = newNodes(j2(tmp)); % column-array
ii = [ii;i2]; jj = [jj;j2]; vv = [vv;v2]; % NCTs1
mL2 = sparse(ii,jj,vv);

%-------------------------------------------------------------------------
% Get error L2
%-------------------------------------------------------------------------
vL2 = uh'*mL2*uh;
vL2 = sqrt(vL2);

end