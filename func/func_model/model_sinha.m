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
    fh.name = 'model_sinha';
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
function valF = findDefF(xx,yy,pa,sub)
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
function exSol = findDefExSol(xx,yy,pa,sub)
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