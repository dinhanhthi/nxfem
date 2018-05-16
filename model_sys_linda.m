%% This is a function handles containing all information about about given
%   in the test case proposed by Linda. It's corresponding to file 
%   main_sys_linda.m. This file includes:
%   - function defF : def of RHS 
%   - function defExSol: def of exact solution
%   - function defPhi: def of interface
%   - (function lambda): the choice of lambda's form
%   - description of domain
%   - (boundary condition)


%% =======================================================================
% DOMAIN: [0,1]X[0,1], interface: circle ((0.5,0.5);r0) inside domain
%-------------------------------------------------------------------------
% MODELS:
% Note that we find tw instead of w (w=beta*tw)
% r=sqrt((x-.5)^2+(y-.5)^2)
%-------------------------------------------------------------------------
% Equation of tw:
% -nabla(beta*nabla(tw)) = -16r^2   in Omega
% [tw]=[beta*nabla_n(tw)]=0     on Gamma
% tw=twex   on \partial\Omega
%-------------------------------------------------------------------------
% Equation of u:
% -nabla(alpha*nabla(u)) + beta*v*g(u) = -4+(wex-alpha*uex)*g(uex) in Omega
% [u]=[alpha*nabla_n(u)]=0 on Gamma
% u=uex on \partial\Omega
%-------------------------------------------------------------------------
% Equation of v: v = 1/beta*(w-alpha*u)
%-------------------------------------------------------------------------
% EXACT SOLUTION:
% twex = r^4/beta1  if r<=r0
%        r^4/beta2 - r0^4/beta2 + r0^4/beta1  if r>r0
% uex = r^2/alpha1  if r<=r0
%       r^2/alpha2 - r0^2/alpha2 + r0^2/alpha1  if r>r0
%=========================================================================



%% =======================================================================
% SETTING UP
%=========================================================================
function fh = model_sys_linda
    fh.defFtw = @findDefFtw; % RHS of equation of tw
    fh.defFu = @findDefFu; % RHS of equation of u
    fh.defWex = @findDefWex; % wex
    fh.defTWex = @findDefTWex; % twex
    fh.defUex = @findDefUex; % uex
    fh.defVex = @findDefVex; % vex
    fh.defPhi = @findDefPhi; % interface
    fh.domain = @findDomain; % domain setting
    fh.bcTW = @findTypeBCtw; % BCs for tw's equation
    fh.bcU = @findTypeBCu; % BCs for u's equation
    fh.kapTW = @findKapTW; % kappa for tw's equation
    fh.kapU = @findKapU; % kappa for u's equation
    fh.lamTW = @findLamTW; % lambda for tw's equation
    fh.lamU = @findLamU; % lambda for u's equation
    fh.lamU = @findLamU; % lambda for u's equation
end



%% =======================================================================
% DOMAIN
%=========================================================================
function GeoDom = findDomain()
    xDomVal = [0 1 1 0]; % x values of points constructing Omega
    yDomVal = [0 0 1 1]; % corresponding y value
    RectDom = [3,4,xDomVal,yDomVal]'; % rectangular domain "3" with "4" sides
    GeoDom = decsg(RectDom);
end



%% =======================================================================
% BCs
%=========================================================================
function typeBC = findTypeBCtw()
    % Output: the number indicating the type of boundary condition
    % find the list of types in file getErrEachStep.m
    typeBC = 2; % 2: u=uex on whole boundary
end
%-------------------------------------------------------------------------
function typeBC = findTypeBCu()
    % Output: the number indicating the type of boundary condition
    % find the list of types in file getErrEachStep.m
    typeBC = 2; % 2: u=uex on whole boundary
end



%% =======================================================================
% F
%=========================================================================
function valF = findDefFtw(xx,yy,sub,pa)
    % Define right hand side function f
    % Input: coordinate of points + indication subdomain
    % Output: value of phi at points
    rr2 = (xx-0.5).^2 + (yy-0.5).^2; % r^2
    valF = -16*rr2; % -16r^2
end
%-------------------------------------------------------------------------
function valF = findDefFu(xx,yy,sub,pa)
    % Define right hand side function fu
    % Input: coordinate of points + indication subdomain
    % Output: value of F at points

    % uex
    uex1 = findDefUex(xx,yy,1,pa);
    uex2 = findDefUex(xx,yy,2,pa);
    % g(u)
    guex1 = defG(uex1,1); 
    guex2 = defG(uex2,1); 
    % wex
    wex1 = findDefWex(xx,yy,1,pa);
    wex2 = findDefWex(xx,yy,2,pa);
    if sub==1 % in Omg1s
        valF = -4 + (wex1-pa.alp1*uex1).*guex1;
    else % in Omg2
        valF = -4 + (wex2-pa.alp2*uex2).*guex2;
    end
end



%% =======================================================================
% EXACT SOLUTIONS
%=========================================================================
function exSol = findDefWex(xx,yy,sub,pa)
    % Describe the exact solution of w
    % Input: coordinate of points + indication subdomain
    % Output: value of exact solution at points
    rr2 = (xx-0.5).^2 + (yy-0.5).^2; % r^2
    if sub==1 % in Omg1s
        exSol = rr2.^2; % r^4
    else % in Omg2
        exSol = rr2.^2-pa.r0^4 + (pa.bet2/pa.bet1)*pa.r0^4;
    end
end
%-------------------------------------------------------------------------
function exSol = findDefTWex(xx,yy,sub,pa)
    % Describe the exact solution of tw
    % Input: coordinate of points + indication subdomain
    % Output: value of exact solution at points
    rr2 = (xx-0.5).^2 + (yy-0.5).^2; % r^2
    if sub==1 % in Omg1s
        exSol = rr2.^2/pa.bet1;
    else % in Omg2
        exSol = (rr2.^2-pa.r0^4)/pa.bet2 + pa.r0^4/pa.bet1;
    end
end
%-------------------------------------------------------------------------
function exSol = findDefUex(xx,yy,sub,pa)
    % Describe the exact solution of u
    % Input: coordinate of points + indication subdomain
    % Output: value of exact solution at points
    rr2 = (xx-0.5).^2 + (yy-0.5).^2; % r^2
    if sub==1 % in Omg1s
        exSol = rr2/pa.alp1;
    else % in Omg2
        exSol = (rr2-pa.r0^2)/pa.alp2 + pa.r0^2/pa.alp1;
    end
end
%-------------------------------------------------------------------------
function exSol = findDefVex(xx,yy,sub,pa)
    % Describe the exact solution of v
    % Input: coordinate of points + indication subdomain
    % Output: value of exact solution at points
    % uex
    uex1 = findDefUex(xx,yy,1,pa);
    uex2 = findDefUex(xx,yy,2,pa);
    % wex
    wex1 = findDefWex(xx,yy,1,pa);
    wex2 = findDefWex(xx,yy,2,pa);
    if sub==1 % in Omg1s
        exSol = 1/pa.bet1*(wex1 - pa.alp1*uex1);
    else % in Omg2
        exSol = 1/pa.bet2*(wex2 - pa.alp2*uex2);
    end
end



%% =======================================================================
% INTERFACE
%=========================================================================
function valPhi=findDefPhi(xx,yy,pa)
    % Define level set function phi
    % Input: coordinate of points
    % Output: value of phi at points
    valPhi = sqrt((xx-0.5).^2+(yy-0.5).^2)-pa.r0;
end



%% =======================================================================
% KAPPA
%=========================================================================
function kap = findKapTW(areaChildCTs,pa)
    % for seeking w
    % kapi : 1 x nCTs
    nCTs = size(areaChildCTs,2);
    kap.kap1 = zeros(1,nCTs) + pa.bet2/(pa.bet1+pa.bet2);
    kap.kap2 = zeros(1,nCTs) + pa.bet1/(pa.bet1+pa.bet2);
end
%-------------------------------------------------------------------------
function kap = findKapU(areaChildCTs,pa)
    % for seeking u
    % kapi : 1 x nCTs
    nCTs = size(areaChildCTs,2);
    kap.kap1 = zeros(1,nCTs) + pa.alp2/(pa.alp1+pa.alp2);
    kap.kap2 = zeros(1,nCTs) + pa.alp1/(pa.alp1+pa.alp2);
end


%% =======================================================================
% LAMBDA (penalty term)
%=========================================================================
function lam = findLamTW(areaChildCTs,pa)
    % for seeking w
    % lam: 1 x nCTs
    nCTs = size(areaChildCTs,2);
    lam = zeros(1,nCTs) + 4*pa.lamHw*pa.bet1*pa.bet2/(pa.bet1+pa.bet2);
end
%-------------------------------------------------------------------------
function lam = findLamU(areaChildCTs,pa)
    % for seeking u
    % lam: 1 x nCTs
    nCTs = size(areaChildCTs,2);
    lam = zeros(1,nCTs) + 4*pa.lamHu*pa.alp1*pa.alp2/(pa.alp1+pa.alp2);
end