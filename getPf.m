function P = getPf(tris,CT,msh,pa,defF)
% get all necessary P for get load vector
% case: int_Omg f(x,y)*phi
% related file: getLf.m
% Input:
% Ouput: - P.NC1, P.NC2: for NCTs
%        - P.CTW1, P.CTW2, P.CTP1, P.CTP2: for CTs

NCTs1=tris.NCTs1; NCTs2=tris.NCTs2; CTs=tris.CTs;
iPs=CT.iPs; typeCTs=CT.type; nodeCTs=CT.nodes;

%% For the whole triangle casse
    function P = getPWhole(Ts,defF,pa,msh,sub)
    % Get P in int_Omg P*phi for both NCTs, CTs (whole triangle case)
    % This P: f(x,y)
    % Related files: function defF in model_sys_linda.m
    % Input: 
    % Output: matrix P: nTs x nwt
    nTs = size(Ts,2);
    points = msh.p;
    dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
    [wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
    nwt = size(wt,2); % number of Gaussian points
    P = zeros(nTs,nwt);
    for t=1:nTs
        triangle = Ts(:,t);
        v1 = points(:,triangle(1)); % vertex 1
        v2 = points(:,triangle(2)); % vertex 2
        v3 = points(:,triangle(3)); % vertex 3
        for k=1:nwt
            [xk,yk] = getCoorSTD(pt(:,k),v1,v2,v3);
            P(t,k) = defF(xk,yk,sub,pa);
        end
    end
    end

%% For the part triangle case
    function P = getPPart(CTs,iPs,typeCTs,nodeCTs,defF,pa,msh,sub)
    % Get P in int_Omg P*phi for both NCTs, CTs (part triangle case)
    % This P: f(x,y)
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
            [xk,yk] = getCoorSTD([xHk,yHk],v1,v2,v3);
            P(t,k) = defF(xk,yk,sub,pa);
        end
    end
    end

%% Get P for NCTs
P.NC1 = getPWhole(NCTs1,defF,pa,msh,1);
P.NC2 = getPWhole(NCTs2,defF,pa,msh,2);

%% Get P for CTs
P.CTW1 = getPWhole(CTs,defF,pa,msh,1);
P.CTW2 = getPWhole(CTs,defF,pa,msh,2);
P.CTP1 = getPPart(CTs,iPs,typeCTs,nodeCTs,defF,pa,msh,1);
P.CTP2 = getPPart(CTs,iPs,typeCTs,nodeCTs,defF,pa,msh,2);

end