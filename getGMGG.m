function A = getGMGG(tris,phi,CT,msh,pa,cp)
% Find the global stiffness matrix for tw (there is only term grad*grad)
% Related file: main_sys_linda, ...
% Status: checked with the old Global matrix function
% This function built from the shape of bilinear form
% Input: - triangles on not cut triangles
%        - information about cut triangles
%	 	 - nodesOmg2GamCTs : to get new nodes in CTs2
% Output: global stiffness matrix Atw

CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
aChild=CT.areaChild; iPs=CT.iPs; uN=CT.uN;
newNodes = msh.newNodes; % convert i to k(i)
kk1 = cp.kk1; kk2 = cp.kk2; % diff coef


%% =======================================================================
% GET TRIPLETS
%=========================================================================

%-------------------------------------------------------------------------
% Term GRAD*GRAD
%-------------------------------------------------------------------------
[i1,j1,v1] = getTripleGGNCTs(NCTs1,kk1,msh); % NCTs1
[i2,j2,v2] = getTripleGGNCTs(NCTs2,kk2,msh); % NCTs2
[iGG,jGG,vGG1,vGG2] = getTripleGGCTs(CTs,aChild,kk1,kk2,msh); % CTs

%-------------------------------------------------------------------------
% 3 terms on interface
%-------------------------------------------------------------------------
% 2 terms grad_n*Phi (sign: "-" for both)
[iGP,jGP,vGP1,vGP2,vGP3,vGP4] = getTripleGPoG(CTs,iPs,uN,msh,pa,cp); 
L = repmat(cp.lambda,4,1);
[iPP,jPP,vPP1,vPP2,vPP3,vPP4] = getTriplePPoG(CTs,iPs,msh,pa,L); 
        % term L*phi*phi (sign +)

%-------------------------------------------------------------------------
% Put into cut triangles cases
%-------------------------------------------------------------------------
it1 = [iGG;iGP;iPP]; jt1 = [jGG;jGP;jPP]; vt1 = [vGG1;vGP1;vPP1]; % A_ij
it2 = it1; jt2 = jt1; vt2 = [vGG2;vGP2;vPP2];  % A_k(i)k(j)
it3 = [iGP;iPP]; jt3 = [jGP;jPP]; vt3 = [vGP3;vPP3]; % A_k(i)j
it4 = it3; jt4 = jt3; vt4 = [vGP4;vPP4]; % A_ik(j)



%% =======================================================================
% GET iG,jG,vG
%=========================================================================

%-------------------------------------------------------------------------
% NCTs1 & CTs1
%-------------------------------------------------------------------------
iG = [i1;it1]; jG = [j1;jt1]; vG = [v1;vt1];

%-------------------------------------------------------------------------
% NCTs2
%-------------------------------------------------------------------------
tmp = ismember(i2,msh.node.CT.omg2); % column-array
i2(tmp) = newNodes(i2(tmp)); % column-array
tmp = ismember(j2,msh.node.CT.omg2); % column-array
j2(tmp) = newNodes(j2(tmp)); % column-array
iG = [iG;i2]; jG = [jG;j2]; vG = [vG;v2];

%-------------------------------------------------------------------------
% CTs2, A_k(i)k(j)
%-------------------------------------------------------------------------
it2 = newNodes(it2); % column-array
jt2 = newNodes(jt2); % column-array
iG = [iG;it2]; jG = [jG;jt2]; vG = [vG;vt2];

%-------------------------------------------------------------------------
% CTs3, A_{k(i)j}
%-------------------------------------------------------------------------
it3 = newNodes(it3); % column-array
iG = [iG;it3]; jG = [jG;jt3]; vG = [vG;vt3];

%-------------------------------------------------------------------------
% CTs4, A_{ik(j)}
%-------------------------------------------------------------------------
jt4 = newNodes(jt4); % column-array
iG = [iG;it4]; jG = [jG;jt4]; vG = [vG;vt4];



%% =======================================================================
% Ghost penalty terms
%=========================================================================
if pa.useGP
    [iGP,jGP,vGP] = getGhostPenalty(CTs,phi,msh,pa,cp);
    iG = [iG;iGP]; jG = [jG;jGP]; vG = [vG;vGP];
end



%% =======================================================================
% GLOBAL MATRIX
%=========================================================================
A = sparse(iG,jG,vG);

end