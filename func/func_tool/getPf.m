function P = getPf(msh,pa,tris,CT,sol,func)
% Get all necessary P for BOTH getL and getGM
% P has form f(w(x,y))*g(u(x,y))*h(x,y)
% Related file: getL* AND getGM*
% Input: - all triangles tris
%        - cut triangle's info
%        - sol.u, sol.w : two variables, 4 components
%        - func.h for h(x,y,pa,sub) function handle
%          func.gu for g(u) function handle
%          func.fw for f(w) function handle
% Ouput: - P.NC1, P.NC2: for NCTs
%        - P.CTW1, P.CTW2, P.CTP1, P.CTP2: for CTs
%       (all are nTs x nwt matrix)


%% Prequirement
NCTs1=tris.NCTs1; NCTs2=tris.NCTs2; CTs=tris.CTs;
% default values for u,w
SolDef.omg1 = ones(msh.nStd,1); SolDef.omg2 = ones(msh.nStd,1);
SolDef.ct1 = ones(msh.nStd,1); SolDef.ct2 = ones(msh.nStd,1);


%% Default functions
if ~isfield(sol,'u')
    uold = SolDef;
    if ~isfield(func,'gu')
       func.gu = @findDefG; 
    end
else
    uold = sol.u;
end
if ~isfield(sol,'w')
    wold = SolDef;
    if ~isfield(func,'fw')
       func.fw = @findDefG; 
    end
else
    wold = sol.w;
end
if ~isfield(func,'h')
    func.h = @findDefH; 
end


%% For the whole triangle casse
    function P = getPWhole(msh,pa,Ts,uold,defGu,wold,defFw,defH,sub)
    % Get P in int_NCTs P*phi for the case of whole triangle
    % This P: f(w(x,y))*g(u(x,y))*h(x,y)
    % Related files: function defF in model_sys_linda.m, defGu.m
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
        uoldT(1:3) = uold(Ts(1:3,t));
        woldT(1:3) = wold(Ts(1:3,t));
        for k=1:nwt
            [shFu,~,~] = getP1shapes(pt(1,k),pt(2,k));
            vuold = uoldT(1)*shFu(1)+uoldT(2)*shFu(2)+uoldT(3)*shFu(3);
            vwold = woldT(1)*shFu(1)+woldT(2)*shFu(2)+woldT(3)*shFu(3);
            [xk,yk] = getCoorSTD(pt(:,k),v1,v2,v3);
            P(t,k) = defGu(vuold)*defFw(vwold)*defH(xk,yk,pa,sub);
        end
    end
    end


%% For the part triangle case
    function P = getPPart(msh,pa,Ts,CT,uold,defGu,wold,defFw,defH,sub)
    % Get P in int_NCTs P*phi for the case of part of triangle
    % This P: f(w(x,y))*g(u(x,y))*h(x,y)
    % Related files: function defF in model_sys_linda.m, defGu.m
    % Input: 
    % Output: matrix P: nTs x nwt
    
    nCTs = size(Ts,2);
    iPs=CT.iPs; typeCTs=CT.type; nodeCTs=CT.nodes;
    points = msh.p;
    dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
    [wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
    nwt = size(wt,2); % number of Gaussian points
    P = zeros(nCTs,nwt);
    for t=1:nCTs
        triangle = Ts(:,t);
        v1 = points(:,triangle(1)); % vertex 1
        v2 = points(:,triangle(2)); % vertex 2
        v3 = points(:,triangle(3)); % vertex 3
        uoldT(1:3) = uold(Ts(1:3,t));
        woldT(1:3) = wold(Ts(1:3,t));
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
            [shFu,~,~] = getP1shapes(xHk,yHk);
            vuold = uoldT(1)*shFu(1)+uoldT(2)*shFu(2)+uoldT(3)*shFu(3);
            vwold = woldT(1)*shFu(1)+woldT(2)*shFu(2)+woldT(3)*shFu(3);
            P(t,k) = defGu(vuold)*defFw(vwold)*defH(xk,yk,pa,sub);
        end
    end
    end


%% Get P for NCTs
P.NC1 = getPWhole(msh,pa,NCTs1,uold.omg1,func.gu,wold.omg1,func.fw,func.h,1);
P.NC2 = getPWhole(msh,pa,NCTs2,uold.omg2,func.gu,wold.omg2,func.fw,func.h,2);


%% Get P for CTs
P.CTW1 = getPWhole(msh,pa,CTs,uold.ct1,func.gu,wold.ct1,func.fw,func.h,1);
P.CTW2 = getPWhole(msh,pa,CTs,uold.ct2,func.gu,wold.ct2,func.fw,func.h,2);
P.CTP1 = getPPart(msh,pa,CTs,CT,uold.ct1,func.gu,wold.ct1,func.fw,func.h,1);
P.CTP2 = getPPart(msh,pa,CTs,CT,uold.ct2,func.gu,wold.ct2,func.fw,func.h,2);



%% Default function handles
function val = findDefH(xx,yy,pa,sub)
    val= 1;
end

function val = findDefG(uu)
    val= 1;
end


end