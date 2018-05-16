function K = getKug(tris,CT,uold,msh,pa,cp,typeG)
% get all necessary K for global matrix
% case: int_Omg K*phi*phi
% This is K for: kk*g(uold)
% related file: getGMgPP.m
% Input:
% Ouput: - K.NC1, K.NC2: for NCTs
%        - K.CTW1, K.CTW2, K.CTP1, K.CTP2: for CTs

NCTs1=tris.NCTs1; NCTs2=tris.NCTs2; CTs=tris.CTs;
iPs=CT.iPs; typeCTs=CT.type; nodeCTs=CT.nodes;
kk1 = cp.kk1; kk2 = cp.kk2; % diff coef

%% For the whole triangle casse
    function K = getKugWhole(Ts,pa,uold,kk,typeG)
    % Get K in int_Omg K*phi*phi for both NCTs, CTs (whole triangle case)
    % This is K for: kk*g(uold)
    % Related files: defG, getTriplePPCTs, getTriplePPNCTs
    % Input: 
    % Output: matrix K: nTs x nwt
    dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
    [wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
    nwt = size(wt,2); % number of Gaussian points
    nTs = size(Ts,2);
    K = zeros(nTs,nwt);
    for t=1:nTs
        uoldT(1:3) = uold(Ts(1:3,t));
        for k=1:nwt
            [shFu,~,~] = getP1shapes(pt(1,k),pt(2,k)); % N_i at quadrature points
            vuold = uoldT(1)*shFu(1)+uoldT(2)*shFu(2)+uoldT(3)*shFu(3);
            guold = defG(vuold,typeG); % g(u), cf. defG.m
            K(t,k) = kk*guold;
        end
    end
    end

%% For the part triangle case
    function K = getKugPart(CTs,iPs,typeCTs,nodeCTs,pa,msh,uold,kk,typeG)
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

%% Get P for NCTs
K.NC1 = getKugWhole(NCTs1,pa,uold.omg1,kk1,typeG); % NCTs1
K.NC2 = getKugWhole(NCTs2,pa,uold.omg2,kk2,typeG); % NCTs2

%% Get P for CTs
K.CTW1 = getKugWhole(CTs,pa,uold.ct1,kk1,typeG); % CTs whole 1
K.CTW2 = getKugWhole(CTs,pa,uold.ct2,kk2,typeG); % % CTs whole 2
K.CTP1 = getKugPart(CTs,iPs,typeCTs,nodeCTs,pa,msh,uold.ct1,kk1,typeG); % CTs part 1
K.CTP2 = getKugPart(CTs,iPs,typeCTs,nodeCTs,pa,msh,uold.ct2,kk2,typeG); % CTs part 2 

end