function vl = getNormAveV(vh,CTs,iPs,msh,cpV,pa)
% Find |{v}|_L2(\Gam)
% Input: - vh & cpV must the same of var
%        - cpV <-- we need kappa_i and kk_i
%        - cut triangles CTs
%        - intersection points iPs
% Output: - a matrix

L = getLAveV(cpV);
[iPP,jPP,vPP1,vPP2,vPP3,vPP4] = getTriplePPoG(CTs,iPs,msh,pa,L);
% IMPORTANT: because we wanna use getTriplePPoG (phi*phi) but there are
% signs minus (-) in lines 3 and 4 whereas in this norm, all signs are +.
% Thus we need to change sign lines 3, 4 of L before applying to this
% function.
newNodes = msh.newNodes;

%-------------------------------------------------------------------------
% Put into cut triangles cases
%-------------------------------------------------------------------------
it1 = iPP; jt1 = jPP; vt1 = vPP1; % A_ij
it2 = it1; jt2 = jt1; vt2 = vPP2;  % A_k(i)k(j)
it3 = iPP; jt3 = jPP; vt3 = vPP3; % A_k(i)j
it4 = it3; jt4 = jt3; vt4 = vPP4; % A_ik(j)

%-------------------------------------------------------------------------
% CTs1
%-------------------------------------------------------------------------
iG = it1; jG = jt1; vG = vt1;

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

mA = sparse(iG,jG,vG);

vl = vh'*mA*vh;
vl = sqrt(vl);
end