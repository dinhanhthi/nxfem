%% This function is just for testing the model given in the thesis of 
% Boiveau (cf. 5.5)
% Related file: main.m, main_eachStep.m
%  This file includes
%   - function defF : def of RHS 
%   - function defExSol: def of exact solution
%   - function defPhi: def of interface
%   - (function lambda): the choice of lambda's form
%   - description of domain
%   - (boundary condition)


%% Model given in article1
% BC: using u=uex
% Test case given in Barrau's thesis (page 113)
% Omg=[0,1]x[0,1], Gam=sqrt((x-1/2)^2+(y-1/2)^2)
%   -nabla\cdot(k\nabla w) = fw in Omg_i
%   [w]=[k\nabla_n w]=0 on Gam
%    w = wex on \part\Omg
% Exact solution :  given in defExSol


%% main function
function fh = model_boiveau
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
xDomVal = [0 1 1 0]; % x values of points constructing Omega
yDomVal = [0 0 1 1]; % corresponding y value
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
if sub==1 % in Omg1s
   valF = -8*pa.kk1*(2*xx.^2-2*xx +2*yy.^2-2*yy +1);
else % in Omg2
   valF = -8*pa.kk2*(2*xx.^2-2*xx +2*yy.^2-2*yy +1);
end
end

%% exact solution
function exSol = findDefExSol(xx,yy,sub,pa)
% Describe the exact solution
% Input: coordinate of points + indication subdomain
% Output: value of exact solution at points

exSol = ((xx-0.5).^2 + (yy-0.5).^2).^2;

end

%% interface
function valPhi=findDefPhi(xx,yy,pa)
% Define level set function phi
% Input: coordinate of points
% Output: value of phi at points
valPhi = sqrt((xx-0.5).^2 + (yy-0.5).^2) - pa.xi;
end


%% =======================================================================
% KAPPA
%=========================================================================
function kap = findKap(areaChildCTs,pa)
    % kapi : 1 x nCTs
    nCTs = size(areaChildCTs,2);
%     kap.kap1 = zeros(1,nCTs) + pa.kk2/(pa.kk1+pa.kk2);
%     kap.kap2 = zeros(1,nCTs) + pa.kk1/(pa.kk1+pa.kk2);
    kap.kap1 = zeros(1,nCTs) + 0.5;
    kap.kap2 = zeros(1,nCTs) + 0.5;
end



%% =======================================================================
% LAMBDA (penalty term)
%=========================================================================
% function lam = findLam(areaChildCTs,pa)
%     % lam: 1 x nCTs
%     nCTs = size(areaChildCTs,2);
%     lam = zeros(1,nCTs) + 4*pa.lamH*pa.kk1*pa.kk2/(pa.kk1+pa.kk2);
% end

function lam = findLam(pa,hT)
    % lam: 1 x nCTs
    coef = 4*pa.lamH*pa.kk1*pa.kk2/(pa.kk1+pa.kk2);
    lam = coef./hT;
end