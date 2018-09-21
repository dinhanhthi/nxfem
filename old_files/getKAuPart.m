function K = getKAuPart(CTs,iPs,typeCTs,nodeCTs,pa,msh,uold,kk,typeG)
% get K in int_CTs(K*phi*phi) for Part triangle
% This is K for: kk*g(uold)
% Related files: defG, getTriplePPCTs
% Input: - triangles (almost CTs)
%        - uold : nstd x 1
%        - kk: (almost) diffusion coefficients
%        - typeG: type of function g(u)
% Output: matrix K (nCTs x nwt)

points = msh.p;
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points
nCTs = size(CTs,2);
K = zeros(nCTs,nwt);
for t=1:nCTs
    uoldT(1:3) = uold(CTs(1:3,t));
    triangle = CTs(:,t); % info of triangle t-th
    v1 = points(:,triangle(1)); % vertex 1
    v2 = points(:,triangle(2)); % vertex 2
    v3 = points(:,triangle(3)); % vertex 3
    iP1 = iPs(:,1,t); % 1st intersection point
    iP2 = iPs(:,2,t); % 2nd intersection point
    if typeCTs(t)==2
        rV = points(:,nodeCTs.eachOmg2(1,t)); 
    else
        rV = points(:,nodeCTs.eachOmg1(1,t)); 
    end
    [xiP1h,yiP1h] = getCoorRef(iP1,v1,v2,v3); % iP1 in ref coor
    [xiP2h,yiP2h] = getCoorRef(iP2,v1,v2,v3); % iP2 in ref coor
    [xRvh,yRvh] = getCoorRef(rV,v1,v2,v3); % remaining vertex in ref coor
    for k=1:nwt
        [xHk,yHk] = getCoorSTD(pt(:,k),[xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);
        [shFu,~,~] = getP1shapes(xHk,yHk);
        vuold = uoldT(1)*shFu(1)+uoldT(2)*shFu(2)+uoldT(3)*shFu(3);
        guold = defG(vuold,typeG); % g(u), cf. defG.m
        K(t,k) = kk*guold;
    end
end
end