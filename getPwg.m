function P = getPwg(tris,CT,uold,wS,msh,pa,typeG)
% get all necessary P for get load vector
% case: int_Omg wg(u)*phi
% related file: getLfwg.m
% Input:
% Ouput: - P.NC1, P.NC2: for NCTs
%        - P.CTW1, P.CTW2, P.CTP1, P.CTP2: for CTs

NCTs1=tris.NCTs1; NCTs2=tris.NCTs2; CTs=tris.CTs;
iPs=CT.iPs; typeCTs=CT.type; nodeCTs=CT.nodes;

%% For the whole triangle casse
    function P = getPWhole(Ts,uold,wS,pa,typeG)
    % Get P in int_Omg P*phi for both NCTs, CTs (whole triangle case)
    % This P: w*g(u)
    % Related files: function defF in model_sys_linda.m
    % Input: 
    % Output: matrix P: nTs x nwt
    nTs = size(Ts,2);
    dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
    [wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
    nwt = size(wt,2); % number of Gaussian points
    P = zeros(nTs,nwt);
    for t=1:nTs
        uoldT(1:3) = uold(Ts(1:3,t));
        wST(1:3) = wS(Ts(1:3,t));
        for k=1:nwt
            [shFu,~,~] = getP1shapes(pt(1,k),pt(2,k));
            vuold = uoldT(1)*shFu(1)+uoldT(2)*shFu(2)+uoldT(3)*shFu(3);
            guold = defG(vuold,typeG); % g(uold)
            vwS = wST(1)*shFu(1)+wST(2)*shFu(2)+wST(3)*shFu(3); % w
            P(t,k) = vwS*guold;
        end
    end
    end

%% For the part triangle case
    function P = getPPart(CTs,iPs,typeCTs,nodeCTs,uold,wS,pa,msh,typeG)
    % Get P in int_Omg P*phi for both NCTs, CTs (part triangle case)
    % This P: w*g(u)
    % Related files: function defF in model_sys_linda.m
    % Input: 
    % Output: matrix P: nTs x nwt
    nCTs = size(CTs,2);
    points = msh.p;
    dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
    [wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
    nwt = size(wt,2); % number of Gaussian points
    P = zeros(nCTs,nwt);
    for t=1:nCTs
        triangle = CTs(:,t);
        uoldT = uold(CTs(1:3,t));
        wST(1:3) = wS(CTs(1:3,t));
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
        [xiP1h,yiP1h] = getCoorRef(iP1,v1,v2,v3); % cut point 1 in ref coor
        [xiP2h,yiP2h] = getCoorRef(iP2,v1,v2,v3); % cut point 2 in ref coor
        [xRvh,yRvh] = getCoorRef(rV,v1,v2,v3); % remaining vertex in ref coor
        for k=1:nwt
            [xHk,yHk] = getCoorSTD(pt(:,k),[xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);
            [shFu,~,~] = getP1shapes(xHk,yHk);
            vuold = uoldT(1)*shFu(1)+uoldT(2)*shFu(2)+uoldT(3)*shFu(3);
            guold = defG(vuold,typeG); % g(uold)
            vwS = wST(1)*shFu(1)+wST(2)*shFu(2)+wST(3)*shFu(3); % w
            P(t,k) = vwS*guold;
        end
    end
    end

%% Get P for NCTs
P.NC1 = getPWhole(NCTs1,uold.omg1,wS.omg1,pa,typeG);
P.NC2 = getPWhole(NCTs2,uold.omg2,wS.omg2,pa,typeG);

%% Get P for CTs
P.CTW1 = getPWhole(CTs,uold.ct1,wS.ct1,pa,typeG);
P.CTW2 = getPWhole(CTs,uold.ct2,wS.ct2,pa,typeG);
P.CTP1 = getPPart(CTs,iPs,typeCTs,nodeCTs,uold.ct1,wS.ct1,pa,msh,typeG);
P.CTP2 = getPPart(CTs,iPs,typeCTs,nodeCTs,uold.ct2,wS.ct2,pa,msh,typeG);

end