%% This is a function handles containing all information about Sinha's
% related files: main.m, main_eachStep.m
%   model. This file includes
%   - function defF : def of RHS 
%   - function defExSol: def of exact solution
%   - function defPhi: def of interface
%   - (function lambda): the choice of lambda's form
%   - description of domain
%   - (boundary condition)

%% Sinha's model
% Test case given in Sinha's paper (STRAIGHT LINE)
% Omg=[0,2]x[0,1], Gam={1}x[0,1]
%   -nabla\cdot(k\nabla u) = f in Omg_i
%   [u]=[k\nabla_n u]=0 on Gam
%   u=0 on \part\Omg
% Exact solution :  given in defExSol.m
% k1 = 2k2 (obligatory)


%% main function
function fh = model_sinha
fh.defF = @findDefF;
fh.defExSol = @findDefExSol;
fh.defPhi = @findDefPhi;
fh.domain = @findDomain;
fh.bc = @findTypeBC;
fh.kap = @findKap; % kappa
fh.lam = @findLam; % lambda 
end


%% =======================================================================
% DOMAIN
%=========================================================================
function GeoDom = findDomain()
% Output: a matrix containing all information of domain
xDomVal = [0 2 2 0]; % x values of points constructing Omega
yDomVal = [0 0 1 1]; % corresponding y value
RectDom = [3,4,xDomVal,yDomVal]'; % rectangular domain "3" with "4" sides
GeoDom = decsg(RectDom);
end


%% =======================================================================
% BCs
%=========================================================================
function typeBC = findTypeBC()
% Output: the number indicating the type of boundary condition
% find the list of types in file getErrEachStep.m
% 1 for u=0, 2 for u=uex
typeBC = 1; % u=0 on whole boundary
end

%% =======================================================================
% F
%=========================================================================
function valF = findDefF(xx,yy,sub,pa)
% Define right hand side function f
% Input: coordinate of points + indication subdomain
% Output: value of phi at points
if sub==1 % in Omg1s
   valF = 2*pa.kk1*(pi)^2.*sin(pi*xx).*sin(pi*yy);
else % in Omg2
   valF = -5*pa.kk2*(pi)^2.*sin(2*pi*xx).*sin(pi*yy);
end
end

%% =======================================================================
% EXACT SOLUTIONS
%=========================================================================
function exSol = findDefExSol(xx,yy,sub,pa)
% Describe the exact solution
% Input: coordinate of points + indication subdomain
% Output: value of exact solution at points
if sub==1 % in Omg1s
   exSol = sin(pi*xx).*sin(pi*yy);
else % in Omg2
   exSol = -sin(2*pi*xx).*sin(pi*yy);
end
end



%% =======================================================================
% INTERFACE
%=========================================================================
function valPhi=findDefPhi(xx,yy,pa)
% Define level set function phi
% Input: coordinate of points
% Output: value of phi at points
valPhi = xx - pa.xi;
end



%% =======================================================================
% KAPPA
%=========================================================================
function kap = findKap(areaChildCTs,pa)
    % kapi : 1 x nCTs
    nCTs = size(areaChildCTs,2);
    kap.kap1 = zeros(1,nCTs) + pa.kk2/(pa.kk1+pa.kk2);
    kap.kap2 = zeros(1,nCTs) + pa.kk1/(pa.kk1+pa.kk2);
end



%% =======================================================================
% LAMBDA (penalty term)
%=========================================================================
function lam = findLam(pa,hT,CTs)
    % lam: 1 x nCTs
    
    hTCT = hT(CTs(5,:));
    coef = 4*pa.lamH*pa.kk1*pa.kk2/(pa.kk1+pa.kk2);
    lam = coef./hTCT;

%     nCTs = size(CTs,2);
%     lam = zeros(1,nCTs) + 4*pa.lamH*pa.kk1*pa.kk2/(pa.kk1+pa.kk2);
end