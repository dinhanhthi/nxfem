function vl = getNormJump1p2(uh,CTs,iPs,msh,pa)
% Find ||[u]||_1/2 = \sum_K hK^-1 \int_GamK [u][v]
% This term like the term int_Gam [u][v], diff is in this case, there is
% h_K^-1, cf. getGMGG.m
% Input: - cut triangles CTs
%        - intersection points iPs
% Output: - a matrix

if ~isempty(CTs)

    hTCTs = msh.hT(CTs(5,:)); % consider only on cut triangles (1 x nCTs)
    hTm1 = 1./hTCTs; % h^{-1}
    L = repmat(hTm1,4,1); % the same L for 4 cases
    [iPP,jPP,vPP1,vPP2,vPP3,vPP4] = getTriplePPoG(CTs,iPs,msh,pa,L);
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

    vl = uh'*mA*uh;
    vl = sqrt(vl);
else
    vl = 0;
end

end