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
function valF = findDefF(xx,yy,pa,sub)
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
function exSol = findDefExSol(xx,yy,pa,sub)
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