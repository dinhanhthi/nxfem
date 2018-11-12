function Av = getGMvAA(tris,phi,vold,wS,CT,msh,pa,cpV,cpW)
% NEED TO MODIFY LATER, FOLLOW getGMgPP.m AND getGMuNewton.m
% Find the global stiffness matrix for v
% grad*grad - {}[] - []{} + {}{} + [][] - gP
% Related file: main_article1, getLfwg
% Status: 
% This function built from the shape of bilinear form
% Input: - triangles on not cut triangles
%        - information about cut triangles
%	 	 - nodesOmg2GamCTs : to get new nodes in CTs2
% Output: global stiffness matrix Av


    function L = getLv(cpV)
        % Get L in int_Gam L*phi*phi
        % L in this case: lam*{}{} + lam*kap1*kap2*[][]
        % Related file: getGMvAA.m, getTriplePPoG.m
        % Input:
        % Output: matrix L: 4 x nCTs
        %           row 1: Aij
        %           row 2: Akikj
        %           row 3: Akij
        %           row 4: Aikj
        kap1 = cpV.kap1; kap2 = cpV.kap2;
        lambda = cpV.lambda;
        nCTs = size(lambda,2);
        L = zeros(4,nCTs);
        L(1,:) = lambda.*kap1;
        L(2,:) = lambda.*kap2;
        % L(3,:) = 0 from zeros
        % L(4,:) = 0 from zeros
    end


CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2; 
aChild=CT.areaChild; iPs=CT.iPs; uN=CT.uN;
newNodes = msh.newNodes; % convert i to k(i)
cp = cpV;
kk1 = cp.kk1; kk2 = cp.kk2; % diff coef


%% =======================================================================
% GET TRIPLETS
%=========================================================================

%-------------------------------------------------------------------------
% Term GRAD*GRAD
%-------------------------------------------------------------------------
[iGG1,jGG1,vGG1] = getTripleGGNCTs(NCTs1,kk1,msh); % NCTs1
[iGG2,jGG2,vGG2] = getTripleGGNCTs(NCTs2,kk2,msh); % NCTs2
[iGGc,jGGc,vGGc1,vGGc2] = getTripleGGCTs(CTs,aChild,kk1,kk2,msh);%CTs

%-------------------------------------------------------------------------
% 3 terms on interface
%-------------------------------------------------------------------------
% 2 terms grad_n*Phi (sign: "-" for both)
[iGP,jGP,vGP1,vGP2,vGP3,vGP4] = getTripleGPoG(CTs,iPs,uN,msh,pa,cp);
% term L*phi*phi (sign +) for: lam*{}{} + lam*kap1*kap2*[][]
L = getLv(cp); % 4 x nCTs
[iPP,jPP,vPP1,vPP2,vPP3,vPP4] = getTriplePPoG(CTs,iPs,msh,pa,L);

%-------------------------------------------------------------------------
% Term -int_Omg lamSys*g(wS-bet/(alp*lamSys)*vold)*phi*phi
%-------------------------------------------------------------------------
typeG = 1; % g(u) as normal, cf. defG.m
K = getKv(tris,CT,vold,wS,msh,pa,cpW,cpV,typeG);
[igPP1,jgPP1,vgPP1] = getTriplePPNCTs(NCTs1,msh,pa,K.NC1); % NCTs1
[igPP2,jgPP2,vgPP2] = getTriplePPNCTs(NCTs2,msh,pa,K.NC2); % NCTs2
[igPPc,jgPPc,vgPPc1,vgPPc2] = getTriplePPCTs(CTs,CT,msh,pa,K); % CTs
% sign "-" in bilinear form
vgPP1=-vgPP1; vgPP2=-vgPP2; vgPPc1=-vgPPc1; vgPPc2=-vgPPc2;

%-------------------------------------------------------------------------
% These terms act like terms grad grad. So we just need to add them to 
%   the terms grad grad  
%-------------------------------------------------------------------------       
iGG1 = [iGG1;igPP1]; jGG1 = [jGG1;jgPP1]; vGG1 = [vGG1;vgPP1]; 
iGG2 = [iGG2;igPP2]; jGG2 = [jGG2;jgPP2]; vGG2 = [vGG2;vgPP2]; 
iGGc = [iGGc;igPPc]; jGGc = [jGGc;jgPPc]; 
vGGc1 = [vGGc1;vgPPc1]; vGGc2 = [vGGc2;vgPPc2];

%-------------------------------------------------------------------------
% Put into cut triangles cases
%-------------------------------------------------------------------------
it1 = [iGGc;iGP;iPP]; jt1 = [jGGc;jGP;jPP]; vt1 = [vGGc1;vGP1;vPP1]; % A_ij
it2 = it1; jt2 = jt1; vt2 = [vGGc2;vGP2;vPP2];  % A_k(i)k(j)
it3 = [iGP;iPP]; jt3 = [jGP;jPP]; vt3 = [vGP3;vPP3]; % A_k(i)j
it4 = [iGP;iPP]; jt4 = [jGP;jPP]; vt4 = [vGP4;vPP4]; % A_ik(j)



%% =======================================================================
% GET iG,jG,vG
%=========================================================================

%-------------------------------------------------------------------------
% NCTs1 & CTs1
%-------------------------------------------------------------------------
iG = [iGG1;it1]; jG = [jGG1;jt1]; vG = [vGG1;vt1];

%-------------------------------------------------------------------------
% NCTs2
%-------------------------------------------------------------------------
tmp = ismember(iGG2,msh.node.CT.omg2); % column-array
iGG2(tmp) = newNodes(iGG2(tmp)); % column-array
tmp = ismember(jGG2,msh.node.CT.omg2); % column-array
jGG2(tmp) = newNodes(jGG2(tmp)); % column-array
iG = [iG;iGG2]; jG = [jG;jGG2]; vG = [vG;vGG2];

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
Av = sparse(iG,jG,vG);

end
