function K = getKv(tris,CT,vold,wS,msh,pa,cpW,cpV,typeG)
% get all necessary K for global matrix
% case: int_Omg K*phi*phi
% This is K for: lamSys*g(wold-bet/(alp*lamSys)*vold)
% related file: getGMvAA.m
% Input:
% Ouput: - K.NC1, K.NC2: for NCTs
%        - K.CTW1, K.CTW2, K.CTP1, K.CTP2: for CTs

NCTs1=tris.NCTs1; NCTs2=tris.NCTs2; CTs=tris.CTs;
iPs=CT.iPs; typeCTs=CT.type; nodeCTs=CT.nodes;
alp1 = cpW.kk1; alp2 = cpW.kk2; % diff coef of w
bet1 = cpV.kk1; bet2 = cpV.kk2; % diff coef of v

%% For the whole triangle casse
    function K = getKvWhole(Ts,pa,vold,wS,bet,alp,typeG)
    % Get K in int_Omg K*phi*phi for both NCTs, CTs (whole triangle case)
    % This is K for: lamSys*g(wold-bet/(alp*lamSys)*vold)
    % Related files: 
    % Input: 
    % Output: matrix K: nTs x nwt
    dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
    [wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
    nwt = size(wt,2); % number of Gaussian points
    nTs = size(Ts,2);
    K = zeros(nTs,nwt);
    for t=1:nTs
        voldT(1:3) = vold(Ts(1:3,t));
        wST(1:3) = wS(Ts(1:3,t));
        for k=1:nwt
            [shFu,~,~] = getP1shapes(pt(1,k),pt(2,k)); % N_i at quadrature points
            vvold = voldT(1)*shFu(1)+voldT(2)*shFu(2)+voldT(3)*shFu(3); % vold
            vwS = wST(1)*shFu(1)+wST(2)*shFu(2)+wST(3)*shFu(3); % wS
            insideG = vwS - bet/(alp*pa.lamSys)*vvold;
            gvold = defG(insideG,typeG); % g(u), cf. defG.m
            K(t,k) = pa.lamSys*gvold;
        end
    end
    end

%% For the part triangle case
    function K = getKvPart(CTs,iPs,typeCTs,nodeCTs,pa,msh,vold,wS,bet,alp,typeG)
    % get K in int_CTs(K*phi*phi) for Part triangle
    % This is K for: lamSys*g(wold-bet/(alp*lamSys)*vold)
    % Related files:
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
        voldT(1:3) = vold(CTs(1:3,t));
        wST(1:3) = wS(CTs(1:3,t));
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
            vvold = voldT(1)*shFu(1)+voldT(2)*shFu(2)+voldT(3)*shFu(3);
            vwS = wST(1)*shFu(1)+wST(2)*shFu(2)+wST(3)*shFu(3); % wS
            insideG = vwS - bet/(alp*pa.lamSys)*vvold;
            guold = defG(insideG,typeG); % g(u), cf. defG.m
            K(t,k) = pa.lamSys*guold;
        end
    end
    end

%% Get P for NCTs
K.NC1 = getKvWhole(NCTs1,pa,vold.omg1,wS.omg1,bet1,alp1,typeG); % NCTs1
K.NC2 = getKvWhole(NCTs2,pa,vold.omg2,wS.omg2,bet2,alp2,typeG); % NCTs2

%% Get P for CTs
K.CTW1 = getKvWhole(CTs,pa,vold.ct1,wS.ct1,bet1,alp1,typeG); % CTs whole 1
K.CTW2 = getKvWhole(CTs,pa,vold.ct2,wS.ct2,bet2,alp2,typeG); % % CTs whole 2
K.CTP1 = getKvPart(CTs,iPs,typeCTs,nodeCTs,pa,msh,vold.ct1,wS.ct1,bet1,alp1,typeG); % CTs part 1
K.CTP2 = getKvPart(CTs,iPs,typeCTs,nodeCTs,pa,msh,vold.ct2,wS.ct2,bet2,alp2,typeG); % CTs part 2 

end