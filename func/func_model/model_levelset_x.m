%% This is a function handles containing all info about model to test with
%   level set implementation 
% Using domain of Barrau_x: a straight line moves to the left at each time step
% This function includes
%   - function defF : def of RHS 
%   - function defPhi: def of interface
%   - description of domain
%   - function velo: def of velocity


%% Model barrau_x
% Omg = [-1,1]x[-1,1]
% Interface: straight line x=c
% velocity vv = (constant 0)
% Equation: \pattial_t phi + vv\dot \nabla phi = 0


%% main function
function fh = model_levelset_x
    fh.defF = @findDefF;
    fh.defPhi = @findDefPhi;
    fh.domain = @findDomain;
    fh.velo = @findVelo;
    fh.veloGrad = @findVeloGrad; % only for testing gradPhi*gradphi
    fh.pa = @findPara;
end

function pa = findPara()
    pa.Tmax = 1;
    pa.dt = 0.05;
end


%% domain
function GeoDom = findDomain()
% Output: a matrix containing all information of domain
%     xDomVal = [-1 1 1 -1]; % x values of points constructing Omega
%     yDomVal = [-1 -1 1 1]; % corresponding y value
    xDomVal = [0 1 1 0]; % x values of points constructing Omega
    yDomVal = [0 0 1 1];
    RectDom = [3,4,xDomVal,yDomVal]'; % rectangular domain "3" with "4" sides
    GeoDom = decsg(RectDom);
end

%% right hand side
function valF = findDefF(xx,yy,sub,pa)
% Define right hand side function f
% Input: coordinate of points + indication subdomain
% Output: value of phi at points
    valF = 0;
end


%% interface
function valPhi=findDefPhi(xx,yy,pa)
% Define level set function phi
% Input: coordinate of points
% Output: value of phi at points
%     pa.xi = -0.21; % interface
    pa.xi = 0.21;
    valPhi = xx - pa.xi;
end


%% velocity
function val=findVelo(xx,yy,tt,sub)
% Define velocity field
% Input: coordinate of points & time
% Output: [x,y]
    if sub==1 % x coor
        val = zeros(1,size(xx,2)) + 0.8;
    else % y coor
        val = zeros(1,size(xx,2));
    end
end

function val=findVeloGrad(xx,yy,tt)
% Define velocity field
% Input: coordinate of points & time
% Output: [x,y]
    val = 0.8*xx;
end