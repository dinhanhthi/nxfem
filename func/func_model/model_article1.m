%% This is a function handles containing all information about model for article 1. 
% It's corresponding to 
%       - file main_article1.m
%       - Note book no.6
%       - file ffpp-article1.edp (freefem++ fitted meshes)
% This file includes:
%   - function defF : def of RHS 
%   - function defExSol: def of exact solution
%   - function defPhi: def of interface
%   - (function lambda): the choice of lambda's form
%   - description of domain
%   - (boundary condition)


%% =======================================================================
% DOMAIN: [-1,1]X[-1,1], interface: circle ((0,0);r0) inside domain
%-------------------------------------------------------------------------
% MODELS:
% r=sqrt(x^2+y^2)-r0 the interface
%-------------------------------------------------------------------------
% Equation of w:
% -nabla(alpha*nabla(w)) = fw   in Omega
% [w]=[alpha*nabla_n(w)]=0     on Gamma
% w=wex   on \partial\Omega
%-------------------------------------------------------------------------
% Equation of v:
% -nabla(beta*nabla(v)) - lam*v*g(u) = fv in Omega
% v=nabla_n(v)=0 on Gamma  <== SPECIAL HERE!!!
% v=vex on \partial\Omega
%-------------------------------------------------------------------------
% Equation of u: u = w - beta/(alpha*lambda)*v
%-------------------------------------------------------------------------
% EXACT SOLUTION:
% uex = r^2/alp1 if r<=r0
%       (r^2-r0^2)/alp2 + r0^2/alp1 if r>r0
% vex = (r^2-r0^2)^2/bet1 if r<=r0  <== SPECIAL HERE!!!
%       0 if r>r0                   <== SPECIAL HERE!!!
% wex = uex + bet/(lam*alp)*vex
%-------------------------------------------------------------------------
% RHS:
% fu = -4+vex*g(uex) = -4 + (r^2-r0^2)^2/bet1*g(r^2/alp1) if r<=r0
%                      -4 if r>r0
% fv = 8(r0^2-2r^2)-lam*vex*g(uex) if r<=r0
%      0 if r<r0
% fw = fu + 1/lam*fv = -4 + 8/lam*(r0^2-2r^2) if r<=r0
%                      -4 if r>r0
%=========================================================================



%% =======================================================================
% SETTING UP
%=========================================================================
function fh = model_article1
    fh.defFw = @findDefFw; % RHS of equation of w
    fh.defFv = @findDefFv; % RHS of equation of v
    fh.defWex = @findDefWex; % wex
    fh.defUex = @findDefUex; % uex
    fh.defVex = @findDefVex; % vex
    fh.defPhi = @findDefPhi; % interface
    fh.domain = @findDomain; % domain setting
    fh.bcW = @findTypeBCw; % BCs for tw's equation
    fh.bcV = @findTypeBCv; % BCs for u's equation
    fh.kapW = @findKap; % kappa for w's equation
    fh.kapV = @findKap; % kappa for u's equation
    fh.lamW = @findLam; % lambda for w's equation
    fh.lamV = @findLam; % lambda for u's equation
end



%% =======================================================================
% DOMAIN
%=========================================================================
function GeoDom = findDomain()
    xDomVal = [-1 1 1 -1]; % x values of points constructing Omega
    yDomVal = [-1 -1 1 1]; % corresponding y value
    RectDom = [3,4,xDomVal,yDomVal]'; % rectangular domain "3" with "4" sides
    GeoDom = decsg(RectDom);
end



%% =======================================================================
% BCs
%=========================================================================
function typeBC = findTypeBCw()
    % Output: the number indicating the type of boundary condition
    % find the list of types in file getErrEachStep.m
    typeBC = 2; % 2: w=wex on whole boundary
end
%-------------------------------------------------------------------------
function typeBC = findTypeBCv()
    % Output: the number indicating the type of boundary condition
    % find the list of types in file getErrEachStep.m
    typeBC = 2; % 2: v=vex on whole boundary
end



%% =======================================================================
% F
%=========================================================================
function valF = findDefFw(xx,yy,pa,sub)
    % Define right hand side function f
    % Input: coordinate of points + indication subdomain
    % Output: value of phi at points
    rr2 = xx.^2 + yy.^2; % r^2
    if sub==1 % in Omg1
        valF = -4 + (8/pa.lamSys)*(pa.r0^2-2*rr2);
    else % in Omg2
        valF = -4;
    end
end
%-------------------------------------------------------------------------
function valF = findDefFv(xx,yy,pa,sub)
    % Define right hand side function fu
    % Input: coordinate of points + indication subdomain
    % Output: value of F at points
    
    rr2 = xx.^2 + yy.^2; % r^2
    % uex
    uex1 = findDefUex(xx,yy,1,pa);
    % vex
    vex1 = findDefVex(xx,yy,1,pa);
    % g(u)
    guex1 = defG(uex1,1); 
    if sub==1 % in Omg1s
        valF = 8*(pa.r0^2-2*rr2)-pa.lamSys*vex1.*guex1;
    else % in Omg2
        valF = 0;
    end
end



%% =======================================================================
% EXACT SOLUTIONS
%=========================================================================
function exSol = findDefWex(xx,yy,pa,sub)
    % Describe the exact solution of w
    % Input: coordinate of points + indication subdomain
    % Output: value of exact solution at points
    
    % uex
    uex1 = findDefUex(xx,yy,pa,1);
    uex2 = findDefUex(xx,yy,pa,2);
    % vex
    vex1 = findDefVex(xx,yy,pa,1);
    vex2 = findDefVex(xx,yy,pa,2);
    if sub==1 % in Omg1s
        exSol = uex1 + pa.bet1/(pa.lamSys*pa.alp1)*vex1;
    else % in Omg2
        exSol = uex2 + pa.bet2/(pa.lamSys*pa.alp2)*vex2;
    end
end
%-------------------------------------------------------------------------
function exSol = findDefUex(xx,yy,pa,sub)
    % Describe the exact solution of u
    % Input: coordinate of points + indication subdomain
    % Output: value of exact solution at points
    
    rr2 = xx.^2 + yy.^2; % r^2
    if sub==1 % in Omg1s
        exSol = rr2/pa.alp1;
    else % in Omg2
        exSol = (rr2-pa.r0^2)/pa.alp2 + pa.r0^2/pa.alp1;
    end
end
%-------------------------------------------------------------------------
function exSol = findDefVex(xx,yy,pa,sub)
    % Describe the exact solution of v
    % Input: coordinate of points + indication subdomain
    % Output: value of exact solution at points
    
    rr2 = xx.^2 + yy.^2; % r^2
    if sub==1 % in Omg1s
        exSol = (rr2-pa.r0^2).^2/pa.bet1;
    else % in Omg2
        exSol = 0;
    end
end



%% =======================================================================
% INTERFACE
%=========================================================================
function valPhi=findDefPhi(xx,yy,pa)
    % Define level set function phi
    % Input: coordinate of points
    % Output: value of phi at points
    valPhi = sqrt(xx.^2+yy.^2)-pa.r0;
end



%% =======================================================================
% KAPPA
%=========================================================================
function kap = findKap(cp,CT,pa)
    % kapi : 1 x nCTs
    
    if pa.useGP
        nCTs = size(CT.areaChild,2);
        kap.kap1 = zeros(1,nCTs) + cp.kk2/(cp.kk1+cp.kk2);
        kap.kap2 = zeros(1,nCTs) + cp.kk1/(cp.kk1+cp.kk2);
    else
        kap.kap1 = cp.kk1*CT.areaChild(2,:) ./ (cp.kk1*CT.areaChild(2,:) + cp.kk2*CT.areaChild(1,:));
        kap.kap2 = cp.kk2*CT.areaChild(1,:) ./ (cp.kk1*CT.areaChild(2,:) + cp.kk2*CT.areaChild(1,:));
    end
end


%% =======================================================================
% LAMBDA (penalty term)
%=========================================================================
function lam = findLam(cp,hTCTs,CT,pa)
    % lam: 1 x nCTs
    nCTs = size(hTCTs, 2);
    
    if pa.useGP
        coef = 4*cp.lamH*cp.kk1*cp.kk2/(cp.kk1+cp.kk2);
        lam = coef./hTCTs;              % belongs to hT
%         lam = zeros(1,nCTs) + coef;   % not belong to hT
    else
        lenAB = zeros(1,nCTs); % nCTs x 1
        lenAB(1,:) = ( (CT.iPs(1,2,:) - CT.iPs(1,1,:)).^2 + (CT.iPs(2,2,:) - CT.iPs(2,1,:)).^2 ) .^(0.5);
        coef = cp.lamH*cp.kk1*cp.kk2 .* lenAB...
            ./(cp.kk2*CT.areaChild(1,:) + cp.kk1*CT.areaChild(2,:));
        lam = coef./hTCTs;              % belongs to hT
%         lam = zeros(1,nCTs) + coef;   % not belong to hT
    end
end