%% This is a function handles containing all information about model for chopp2007
% It's corresponding to 
%       - file main_chopp2007.m
%       - file choppSimpleKap.edp (freefem++)
% This file includes:
%   - function defF : def of RHS 
%   - function defPhi: def of initial interface
%   - (function lambda): the choice of lambda's form
%   - description of domain
%   - (boundary condition)


%% =======================================================================
% DOMAIN: [0,1]X[0,1], inital interface: half of circle ((1/2,0);r0)
%-------------------------------------------------------------------------
% MODELS:
% r=sqrt(x^2+y^2)-r0 the interface
%-------------------------------------------------------------------------
% Equation of u (substrate):
% -nabla(alpha*nabla(u)) = -mu*u  in Omega
% [u]=[alpha*nabla_n(u)]=0     on Gamma
% u = 1e-5   on \partial\Omega_3
% nabla_n u = 0 elsewhere
%-------------------------------------------------------------------------
% Equation of v (potential):
% -nabla(nabla(v)) = -beta*u in Omega
% v=nabla_n(v)=0 on Gamma 
% v=0 on \partial\Omega_3
% nabla_n v = 0 elsewhere
%-------------------------------------------------------------------------
% PARAMETERS:
% alpha = 120 (Omg1), 150 (Omg2)
% beta = 1e6 (Omg1), 0 (Omg2)
% mu = 3.6e6 (Omg1), 0 (Omg2)
%=========================================================================



%% Setting up
function fh = model_chopp2007
    fh.defG = @findDefGu; % g(u)
    fh.defFu = @findDefFu; % RHS of equation of w
    fh.defFv = @findDefFv; % RHS of equation of v
    fh.defPhi = @findDefPhi; % interface
    fh.domain = @findDomain; % domain setting
    fh.bcU = @findTypeBCu; % BCs for u's equation
    fh.bcV = @findTypeBCv; % BCs for u's equation
    fh.kapU = @findKapU; % kappa for u's equation
    fh.kapV = @findKapV; % kappa for u's equation
    fh.lamU = @findLamU; % lambda for tw's equation
    fh.lamV = @findLamV; % lambda for u's equation
    fh.pa = @findPara; % parameters
end


%% Finding parameters
function pa = findPara()
    pa.Tmax = 1;
    pa.dt = 0.05;
end


%% g(u)
function valG = findDefGu(u)
    valG = 1;
end


%% Domain
function GeoDom = findDomain()
    xDomVal = [0 1 1 0]; % x values of points constructing Omega
    yDomVal = [0 0 1 1]; % corresponding y value
    RectDom = [3,4,xDomVal,yDomVal]'; % rectangular domain "3" with "4" sides
    GeoDom = decsg(RectDom);
end


%% BCs
function typeBC = findTypeBCu()
    % Output: the number indicating the type of boundary condition
    % find the list of types in file getErrEachStep.m
    % 1: w=0 on whole boundary, 2:wex, 3: constant on 3
    typeBC = 3; 
end

function typeBC = findTypeBCv()
    % Output: the number indicating the type of boundary condition
    % find the list of types in file getErrEachStep.m
    % 1: w=0 on whole boundary, 2:wex, 3: constant on 3
    typeBC = 3;
end


%% F
function valF = findDefFu(xx,yy,pa,sub)
    % Define right hand side function f
    % Input: coordinate of points + indication subdomain
    % Output: value of phi at points
    valF = 0;
end

function valF = findDefFv(xx,yy,pa,sub)
    % Define right hand side function fu
    % Input: coordinate of points + indication subdomain
    % Output: value of F at points
    valF = 0;
end


%% Interface
function valPhi=findDefPhi(xx,yy,pa)
    % Define level set function phi
    % Input: coordinate of points
    % Output: value of phi at points
    valPhi = sqrt((xx-0.5).^2+yy.^2)-pa.r0;
end


%% Kappa
function kap = findKapU(areaChildCTs,pa)
    % for seeking w
    % kapi : 1 x nCTs
    nCTs = size(areaChildCTs,2);
    kap.kap1 = zeros(1,nCTs) + pa.alp2/(pa.alp1+pa.alp2);
    kap.kap2 = zeros(1,nCTs) + pa.alp1/(pa.alp1+pa.alp2);
end
%-------------------------------------------------------------------------
function kap = findKapV(areaChildCTs,pa)
    % for seeking u
    % kapi : 1 x nCTs
    nCTs = size(areaChildCTs,2);
    kap.kap1 = zeros(1,nCTs) + 0.5;
    kap.kap2 = zeros(1,nCTs) + 0.5;
end


%% Lambda
function lam = findLamU(cp,hT,CTs,pa)
    % lam: 1 x nCTs
    
    hTCT = hT(CTs(5,:));
    coef = 4*pa.lamHu*cp.kk1*cp.kk2/(cp.kk1+cp.kk2);
    lam = coef./hTCT;

%     nCTs = size(CTs,2);
%     lam = zeros(1,nCTs) + 4*pa.lamHu*cp.kk1*cp.kk2/(cp.kk1+cp.kk2);
end
%-------------------------------------------------------------------------
function lam = findLamV(cp,hT,CTs,pa)
    % lam: 1 x nCTs
    
    % belongs hT
    hTCT = hT(CTs(5,:));
    coef = 4*pa.lamHv*cp.kk1*cp.kk2/(cp.kk1+cp.kk2);
    lam = coef./hTCT;

    % not belong hT
%     nCTs = size(CTs,2);
%     lam = zeros(1,nCTs) + 4*pa.lamHv*cp.kk1*cp.kk2/(cp.kk1+cp.kk2);
end