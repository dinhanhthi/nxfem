function P = getMatPrecond(NCTs1,NCTs2,CTs,pa,areaChildCTs,msh)
% NEED TO MODIFY LATER, FOLLOW getGMgPP.m AND getGMuNewton.m
% Get the preconditioner like the idea at page 1072 on article of Zunino2011
% Input: - matrix A
%        - and other neccessary parameters
% Output: a square sparse matrix P
% NOTE: it's done but not checked, come back later!!!

k1=pa.kk1; k2=pa.kk2;
newNodes = msh.newNodes;
nodesOutsideCTs = setdiff(1:msh.ndof,[msh.node.CT.all,msh.nStd:msh.ndof]);

%% get term gradgrad on NCTs
[i1,j1,v1] = getTripleGGNCTs(NCTs1,k1,msh); % NCTs1
[i2,j2,v2] = getTripleGGNCTs(NCTs2,k2,msh); % NCTs2
mNCTs1 = sparse(i1,j1,v1); % in NCTs1
mNCTs2 = sparse(i2,j2,v2); % in NCTs2
[i1,j1,v1] = find(mNCTs1); [i2,j2,v2] = find(mNCTs2); % on NCTs
iG = i1; jG = j1; vG = v1;
% Replace all nodes in i2 and j2 which are also in nodesCTsInOmg2OnGam with
%   the new values
tmp = ismember(i2,msh.node.CT.omg2); % column-array
i2(tmp) = newNodes(i2(tmp)); % column-array
tmp = ismember(j2,msh.node.CT.omg2); % column-array
j2(tmp) = newNodes(j2(tmp)); % column-array
% Add to global matrix
iG = [iG;i2]; jG = [jG;j2]; vG = [vG;v2];
AoutCTs = sparse(iG,jG,vG); % doesn't contain values at nodes in CTs

%%-----------------------------------
[iGG,jGG,vGG1,vGG2] = getTripleGGCTs(CTs,areaChildCTs,k1,k2,msh); % CTs (term GradGrad)
mCTs1 = sparse(iGG,jGG,vGG1);
mCTs2 = sparse(iGG,jGG,vGG2);
[it1,jt1,vt1] = find(mCTs1); [it2,jt2,vt2] = find(mCTs2); % on CTs
iT=it1; jT=jt1; vT=vt1; 
% Replace all nodes in it2, jt2 with the new values 
% (because these nodes are in the cut triangles region)
it2 = newNodes(it2); % column-array
jt2 = newNodes(jt2); % column-array
% Add to global matrix
iT = [iT;it2]; jT = [jT;jt2]; vT = [vT;vt2];
AinCTs = sparse(iT,jT,vT); % contain only values at nodes in CTs


%%---------------------------
diagAinCTs = diag(diag(AinCTs)); % diagonal matrix of AinCTs
P = diagAinCTs;
P(nodesOutsideCTs,nodesOutsideCTs) = AoutCTs(nodesOutsideCTs,nodesOutsideCTs);


end
