%% This is a function handles containing all information about Barrau's
% related files: main.m, main_eachStep.m
%   model. This file includes
%   - function defF : def of RHS 
%   - function defExSol: def of exact solution
%   - function defPhi: def of interface
%   - (function lambda): the choice of lambda's form
%   - description of domain
%   - (boundary condition)


%% Barrau's model
% Test case given in Barrau's paper (STRAIGHT LINE)
% BC: using u=uex
% Test case given in Barrau's thesis (page 44)
% Omg=[-1,1]x[-1,1], Gam={xi}x[-1,1]
%   -nabla\cdot(k\nabla u) = -2 in Omg_i
%   [u]=[k\nabla_n u]=0 on Gam
%    u = uex on \part\Omg
% Exact solution :  given in defExSol


%% main function
function fh = model_barrau_x
fh.defF = @findDefF;
fh.defExSol = @findDefExSol;
fh.defPhi = @findDefPhi;
fh.domain = @findDomain;
fh.bc = @findTypeBC;
fh.kap = @findKap; % kappa
fh.lam = @findLam; % lambda 
end


%% domain
function GeoDom = findDomain()
% Output: a matrix containing all information of domain
xDomVal = [-1 1 1 -1]; % x values of points constructing Omega
yDomVal = [-1 -1 1 1]; % corresponding y value
RectDom = [3,4,xDomVal,yDomVal]'; % rectangular domain "3" with "4" sides
GeoDom = decsg(RectDom);
end

function typeBC = findTypeBC()
% Output: the number indicating the type of boundary condition
% find the list of types in file getErrEachStep.m
typeBC = 2; % 2: u=uex on whole boundary
end

%% right hand side
function valF = findDefF(xx,yy,sub,pa)
% Define right hand side function f
% Input: coordinate of points + indication subdomain
% Output: value of phi at points
valF = -2;
end

%% exact solution
function exSol = findDefExSol(xx,yy,sub,pa)
% Describe the exact solution
% Input: coordinate of points + indication subdomain
% Output: value of exact solution at points
if sub==1 % in Omg1s
   exSol = xx.^2/pa.kk1;
else % in Omg2
   exSol = (xx.^2-pa.xi^2)/pa.kk2 + pa.xi^2/pa.kk1;
end
end

%% interface
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
function lam = findLam(areaChildCTs,pa)
    % lam: 1 x nCTs
    nCTs = size(areaChildCTs,2);
    lam = zeros(1,nCTs) + 4*pa.lamH*pa.kk1*pa.kk2/(pa.kk1+pa.kk2);
end