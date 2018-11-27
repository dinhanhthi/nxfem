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
% Test case given in Barrau's thesis (CIRCLE INSIDE DOMAIN)
% BC: using u=uex
% Test case given in Barrau's thesis (page 113)
% there is a wrong point in her thesis, it's r^2 instead of r^4
% Omg=[-1,1]x[-1,1], Gam=sqrt(x^2+y^2)
%   -nabla\cdot(k\nabla u) = -4 in Omg_i
%   [u]=[k\nabla_n u]=0 on Gam
%    u = uex on \part\Omg
% Exact solution :  given in defExSol


%% main function
function fh = model_barrau_c
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
function valF = findDefF(xx,yy,pa,sub)
% Define right hand side function f
% Input: coordinate of points + indication subdomain
% Output: value of phi at points
valF = -4;
end

%% exact solution
function exSol = findDefExSol(xx,yy,pa,sub)
% Describe the exact solution
% Input: coordinate of points + indication subdomain
% Output: value of exact solution at points
rr2 = xx.^2+yy.^2; % r^2
if sub==1 % in Omg1s
   exSol = rr2/pa.kk1;
else % in Omg2
   exSol = (rr2-pa.xi^2)/pa.kk2 + pa.xi^2/pa.kk1;
end
end

%% interface
function valPhi=findDefPhi(xx,yy,pa)
% Define level set function phi
% Input: coordinate of points
% Output: value of phi at points
valPhi = sqrt(xx.^2+yy.^2)-pa.xi;
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
        kap.kap1 = cp.kk2*CT.areaChild(1,:) ./ (cp.kk1*CT.areaChild(2,:) + cp.kk2*CT.areaChild(1,:));
        kap.kap2 = cp.kk1*CT.areaChild(2,:) ./ (cp.kk1*CT.areaChild(2,:) + cp.kk2*CT.areaChild(1,:));
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