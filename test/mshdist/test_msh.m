
addpath(genpath('../../func')); % add all necessary functions


%% =======================================================================
% PARAMETERS
% Note that parameters are different for equations of w and u
%=========================================================================

%-------------------------------------------------------------------------
% Fixed parameters
%-------------------------------------------------------------------------
pa.degP1D = 3; % Gaussian quadrature points in 1D (polinomial functions)
pa.degP2D = 4; % Gaussian quadrature points in 2D (polinomial functions)
pa.degN = 8; % Gaussian quadrature points in 2D (non-polynomial functions)
% degree-#OfPoints : 1-1, 2quick -3, 3-4, 4-6, 5-7, 6-12,
%                    7-13, 8-16, 9-19, 10-25, 11-27, 12-33
pa.tol = eps(1e3); % tolerance, 1e-14
model = model_barrau_x; % Becker's test case with interface: x=x_0
% typeG: modify directly in defG.m



%-------------------------------------------------------------------------
% Settings
%-------------------------------------------------------------------------
nSeg = 15; % mesh settings
velo = [.01;0]; % Velocity (grad of potential in other cases)
pa.xi = -0.21; % initial interface
pa.reguMesh = 0; % use regular mesh or not?
pa.smallCut = 0; % ignore small-support basis (1=ignore,0=no)
pa.tH = 100; % to find the small support using (20) or (21) in arnold 2008


%% =======================================================================
% DOMAIN
%=========================================================================
GeoDom = model.domain(); % domain

%-------------------------------------------------------------------------
% Mesh setting up
%-------------------------------------------------------------------------
if ~pa.reguMesh % not regular mesh?
    hEdgeMax = 2/nSeg;
    [points,edges,triangles] = initmesh(GeoDom,'hmax',hEdgeMax); %irregular
else
    [points,edges,triangles] = poimesh(GeoDom,nSeg,nSeg); % regular
end
msh.p = points; msh.t = triangles; msh.e = edges; % save to msh
x = points(1,:); % x-coordinate of points
y = points(2,:); % y-coordinate of points
% diameter (longest side) of each triangle: 1 x nTs
msh.hT = getDiam(msh); % 1 x number of triangles
msh.hTmax = max(msh.hT); % maximum of all diameters
hTmax = msh.hTmax;

%-------------------------------------------------------------------------
% Level set function (INITIAL)
%-------------------------------------------------------------------------
phi = model.defPhi(x,y,pa); % 1 x number of points (row array)
phi(abs(phi)<pa.tol)=0; % find phi which are very small (~0) and set to 0

pathfile = '';

mshdist_w_mesh(msh,pathfile,'phi0'); % mesh
mshdist_w_sol(msh,phi,pathfile,'phi0'); % original phi

phi1 = phi.^3; % modify phi

mshdist_w_mesh(msh,pathfile,'phi1'); % mesh
mshdist_w_sol(msh,phi1,pathfile,'phi1'); % modified phi

mshdist_w_mesh(msh,pathfile,'phi1_old'); % mesh
system('cp phi1.sol phi1_old.sol');
status = system('mshdist phi1');

phinew =  mshdist_r_sol(phi,pathfile,'phi1'); % redistanced phi

% nf = 0; % reset every loop to be sure uh, vh plotted on the same figure
% plotNXFEM(msh,iPs,nf,phi,'withMesh',false,'title','phi ori','dim',2,'export',false); % original phi
% plotNXFEM(msh,iPs,nf,phi,'withMesh',false,'title','phi mod','dim',2,'export',false); % modified phi
% plotNXFEM(msh,iPs,nf,phi,'withMesh',false,'title','phi redist','dim',2,'export',false); % redistanced phi
