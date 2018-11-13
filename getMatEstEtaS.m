function metaS = getMatEstEtaS(iPs,CTs,phi,pa,msh)
% Try to find a suitable value for lambda by using the of posteriori error
%   estimator. The form of estimator is given in page 33, Barrau's thesis.
% This function computes the matrix to find the value of term etaS
% Input: all neccessary inputs
% Output: matrix to be used in finding value of etaS

newNodes = msh.newNodes; % convert i to k(i)

%% ========================================================
% GET TRIPLETS
% all are column-arrays
% =========================================================
% [i1,j1,v1] = getMatEstEtaSNCE(omg1NCTs,1,phi,msh,pa); % not-cut edges in Omg1
% [i2,j2,v2] = getMatEstEtaSNCE(omg2NCTs,2,phi,msh,pa); % not-cut edges in Omg2
[it,jt,vt1,vt2,vt3,vt4] = getMatEstEtaSonGam(CTs,iPs,msh,pa); % interface segments
[ic,jc,vc1,vc2] = getMatEstEtaSCE(CTs,msh,pa,phi); % cut edges


                                    
%% ========================================================
% BUID MATRIX
% =========================================================
% omg1NCTs & CTs1
% -------------------------
% ii = [i1;it]; jj = [j1;jt]; vv = [v1;vt1];
ii = it; jj = jt; vv = vt1; % don't consider omg1NCTs

% CTs2
% -------------------------
itmp = newNodes(it); jtmp = newNodes(jt);
ii = [ii;itmp]; jj = [jj;jtmp]; vv = [vv;vt2]; 

% omg2NCTs
% -------------------------
% tmp = ismember(i2,nodesCTsInOmg2OnGam); % column-array
% i2(tmp) = newNodes(i2(tmp)); % column-array
% tmp = ismember(j2,nodesCTsInOmg2OnGam); % column-array
% j2(tmp) = newNodes(j2(tmp)); % column-array
% ii = [ii;i2]; jj = [jj;j2]; vv = [vv;v2]; % omg2NCTs

% CTs3
% -------------------------
it3 = newNodes(it); % column-array
ii = [ii;it3]; jj = [jj;jt]; vv = [vv;vt3];

% CTs4
% -------------------------
jt4 = newNodes(jt); % column-array
ii = [ii;it]; jj = [jj;jt4]; vv = [vv;vt4];

% cut edges 1
% -------------------------
ii = [ii;ic]; jj = [jj;jc]; vv = [vv;vc1];

% cut edges 2
% -------------------------
ic2 = newNodes(ic); jc2 = newNodes(jc);
ii = [ii;ic2]; jj = [jj;jc2]; vv = [vv;vc2];

metaS = sparse(ii,jj,vv);
end