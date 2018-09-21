%% This is a function handles containing all info about model to test with
%   level set implementation 
% File: implementation standard level set method - niklas johansson.pdf
%   (section 4.1)
% This function includes
%   - function defF : def of RHS 
%   - function defPhi: def of interface
%   - description of domain
%   - function velo: def of velocity



%% Model vortex
% Given in section 4.1 (page 5)
% Omg = [0,1]x[0,1]
% Interface: circle center (0.5;0.75) with radius 0.15
% velocity vv = (u v)
%   u = 2sin(2pi*y)sin^2(pi*x)cos(pi*t)
%   v = -2sin(2pi*x)sin^2(pi*y)cos(pi*t)
% Equation: \pattial_t phi + vv\dot \nabla phi = 0



%% NOTE:
% When t=1, phi comes back to the initial position
% When t=0.5, direction of vertex change 
% time step delt/grid size h = constant
% inside circle phi takes positive and vice versa
% boundary condition: Arnold Book p.221, due to velo=0 on boundary, we don't need BCs for phi



%% main function
function fh = model_levelset_vortex
    fh.name = 'model_level_vortex';
    fh.defF = @findDefF;
    fh.defPhi = @findDefPhi;
    fh.domain = @findDomain;
    fh.velo = @findVelo;
    fh.pa = @findPara; % parameters
end



%% Find parameters
function pa = findPara()
    pa.Tmax = 1;
    pa.dt = 0.01;
%     pa.dt = 0.1;
end



%% domain
function GeoDom = findDomain()
% Output: a matrix containing all information of domain
    xDomVal = [0 1 1 0]; % x values of points constructing Omega
    yDomVal = [0 0 1 1]; % corresponding y value
    RectDom = [3,4,xDomVal,yDomVal]'; % rectangular domain "3" with "4" sides
    GeoDom = decsg(RectDom);
end



%% right hand side
function valF = findDefF(xx,yy,pa,sub)
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
    pa.xi = 0.15;
    valPhi = sqrt((xx-0.5).^2+(yy-0.75).^2)-pa.xi;
end



%% velocity
function val=findVelo(xx,yy,tt,sub)
% Define velocity field
% Input: coordinate of points & time
% Output: [x,y]
    if sub==1 % x coor
        val = 2*sin(2*pi*yy).*(sin(pi*xx)).^2.*cos(pi*tt);
    else % y coor
        val = -2*sin(2*pi*xx).*(sin(pi*yy)).^2.*cos(pi*tt);
    end
end