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
% ------------------------------------------------------------------------
% Last modified: 08-11-2018
% ------------------------------------------------------------------------


%% =======================================================================
% DOMAIN: [0,0.5]X[0,0.5], inital interface: half of circle ((1/4,0);r0)
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
    fh.name = 'model_chopp2007';
    fh.defFu = @findDefFu;      % RHS of equation of u
    fh.defFv = @findDefFv;      % RHS of equation of v
    fh.defPhi = @findDefPhi;    % interface
    fh.domain = @findDomain;    % domain setting
    fh.bcU = @findTypeBCu;      % BCs for u's equation
    fh.bcV = @findTypeBCv;      % BCs for u's equation
%     fh.kapU = @findKapU;        % kappa for u's equation
    fh.kapU = @findKapV;        % kappa for u's equation
    fh.kapV = @findKapV;        % kappa for u's equation
%     fh.lamU = @findLamU;      % lambda for u's equation
    fh.lamU = @findLamV;        % lambda for u's equation
    fh.lamV = @findLamV;        % lambda for u's equation
    fh.pa = @findPara;          % parameters
end



%% Domain
function GeoDom = findDomain()
    xDomVal = [0 0.5 0.5 0];                % x values of points constructing Omega
    yDomVal = [0 0 0.5 0.5];                % corresponding y value
    RectDom = [3,4,xDomVal,yDomVal]';   % rectangular domain "3" with "4" sides
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
% (only constant or given function f(x,y)
% Don't contains something like g(u) where we're computing v
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
    
    if ~pa.phiNew
%     valPhi = sqrt((xx-0.25).^2 + yy.^2) - pa.r0; % semi circle (no need to initialize)
    %     valPhi = pa.a^2*(xx-0.5).^2 + yy.^2/(pa.a^2) - pa.r0^2; % need to use mshdist to initilize to be a signed distance function
    else
        valPhi = (yy-pa.phiHeight) + pa.phiNoise*cos(4*pi*xx/0.5); % need to initialize
    end
end


%% Kappa
%-------------------------------------------------------------------------
function kap = findKapV(cp,CT,pa)
    % for seeking u and v
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


%% Lambda
%-------------------------------------------------------------------------
function lam = findLamV(cp,hTCTs,CT,pa)
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