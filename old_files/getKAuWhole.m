function K = getKAuWhole(NCTs,pa,uold,kk,typeG)
% get K in int_NCTs(K*phi*phi) and K in int_CTs for whole
% This is K for: kk*g(uold)
% Related files: defG, getTriplePPNCTs, getTriplePPCTs
% Input: - triangles (NCTs or CTs)
%        - uold: nStd x 1
%        - diff coeff kk
%        - type of g(u)
% Output: matrix K (nNCTs x nwt)

dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points
nNCTs = size(NCTs,2);
K = zeros(nNCTs,nwt);
for t=1:nNCTs
    uoldT(1:3) = uold(NCTs(1:3,t));
    for k=1:nwt
        [shFu,~,~] = getP1shapes(pt(1,k),pt(2,k)); % N_i at quadrature points
        vuold = uoldT(1)*shFu(1)+uoldT(2)*shFu(2)+uoldT(3)*shFu(3);
        guold = defG(vuold,typeG); % g(u), cf. defG.m
        K(t,k) = kk*guold;
    end
end
end