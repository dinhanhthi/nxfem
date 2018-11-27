%% add path of functions
addpath(genpath('func')); % add all necessary functions



%% Fixed parameters
pa.degP1D = 3; % Gaussian quadrature points in 1D (for polinomial functions)
pa.degP2D = 4; % Gaussian quadrature points in 2D (for polinomial functions)
pa.degN = 8; % Gaussian quadrature points in 2D (for non-polynomial functions)
% degree-#OfPoints : 1-1, 2-3, 3-4, 4-6, 5-7, 6-12, 7-13, 8-16, 9-19, 10-25, 11-27, 12-33
pa.tol = eps(1e3); % tolerance, 1e-14



%% model 
mdl = 3; 
switch mdl
    case 1 % sinha (cf. model_sinha.m)
        model=model_sinha; % file sinha.m
        cp.kk1 = 1; cp.kk2 = 0.5; % diffusion coefficients (kk1 must be = 2kk2)
        pa.xi = 1.5; % DON'T CHANGE!!!
    case 2 % barrau_x (p.44) STRAIGHT LINE [-1,1]x[-1,1]
        model=model_barrau_x; % file barrau_x.m
%         cp.kk1 = 1; cp.kk2 = 1000; % for barrau case
        cp.kk1 = 1; cp.kk2 = 1; % testing
        pa.xi = 0.31; % only used for barrau's model, between (-1,1)
    case 3 % barrau_r (p.47) PART OF CIRCLE [0,1]x[0,1]
        model=model_barrau_r; % file barrau_r.m
        cp.kk1 = 1; cp.kk2 = 100; % for barrau case
        pa.xi = 0.61; % only used for barrau's model, between (0,1)
    case 4 % barrau_c (p.113) CIRCLE INSIDE DOMAIN [-1,1]x[-1,1]
        model=model_barrau_c;
        cp.kk1 = 1; cp.kk2 = 1000;
        pa.xi = 0.71;
    case 5 % only w in article1
        model=model_w;
%         cp.kk1=1; cp.kk2=100;
        cp.kk1=1; cp.kk2=1;
        pa.xi=0.6;
        pa.lamSys=1;
    case 6 % test case in Boiveau's thesis (cf. 5.5)
        % note that, this is for nonsysmmetric + penalty free NXFEM method
        model = model_boiveau;
        cp.kk1=1; cp.kk2 = 100;        
%         cp.kk1=10; cp.kk2 = 1;
        pa.xi=0.3;
end
GeoDom = model.domain(); % domain
pa.kk1 = cp.kk1; pa.kk2 = cp.kk2; % just for defF & defEx



%% SETTINGS
findCR = 0; % wanna find convergence rate? 1 or 0
%     numSegCR = [32, 64, 128, 256]; % only works with findCR=1
%     numSegCR = [32, 62, 94, 128];
    numSegCR = [16,32,64];
    showPlotCR = 0; % show plot of convergence (for findCR=1)
numSegPlot = 5; % only for plotting, findCR=0
savePlot = 0; % 1 = export figures to files (and without plotting them)
    showPlot = 0; % wanna plot or not the solution? (JUST FOR savePlot=0)
    nf = 0; % counter of figures (plot each plot in a separated figure)
pa.smallCut = 0; % ignore small-support basis (1=ignore,0=no)
    pa.tH = 10; % to find the small support using (20) or (21) in arnold 2008
reguMesh = 1; % regular or irregular mesh?
    
cp.lamH = 1e5; % penalty coefficient

% ghost penalty
pa.useGP = 0; % wanna use ghost penalty term?
    pa.gam1 = 1e-5; % parameter for 1st term
    pa.gam2 = 1e-5; % parameter for 2nd term

    
    
%% Dependent parameters
if findCR == 1
    showPlot = 0;
    numSeg = numSegCR;
    disp('Find the convergence rate...');
else
    numSeg = numSegPlot;
    disp('Only ONE step...');
end
nStep = size(numSeg,2);


%% For finding CR
hTarray = zeros(1,nStep);
nDOFsArray = zeros(1,nStep);
errArray = zeros(2,nStep);



%% Foor loops: find the convergence rate if needed
for z = 1:nStep
    
    disp('-----------------------------');
    fprintf('nSeg = %d\n',numSeg(z));
    
    %% Get mesh info
    hMax = 2/numSeg(z);
    if ~reguMesh
        [points,edges,triangles] = initmesh(GeoDom,'hmax',hMax); % irregular
    else
        [points,edges,triangles] = poimesh(GeoDom,numSeg(z),numSeg(z)); % regular
    end
    msh.p = points; msh.t = triangles; msh.e = edges;
    x = points(1,:); % x-coordinate of points
    y = points(2,:); % y-coordinate of points
    msh.hT = getDiam(msh); % 1 x number of triangles  (diameter (longest side) of each triangle)
    msh.hTmax = max(msh.hT); % maximum of all diameters
    
    
    %% Level set function
    defPhi = model.defPhi; % function handle
    phi = defPhi(x,y,pa); % 1 x number of points (row array)
    phi(abs(phi)<pa.tol)=0; % find phi which are very small (~0) and set to 0

    
    %% Get all triangles
    tris = getTriangles(phi,msh,pa);
    CTs=tris.CTs;
    
    
    %% CT's info
    CT = getInfoCTs(CTs,phi,msh,pa);
    nodeCTs=CT.nodes; areaChildCTs=CT.areaChild; iPs=CT.iPs;
    
    
    %% Small cut
    if pa.smallCut
        tic;time=0;
        fprintf('Removing small cut triangles... ')
        [tris,CT] = findSmallPhi_after(msh,pa,phi,tris,CT);
%         clear CTs NCTs NCTs2 nodeCTs areaChildCTs iPs; % just in case
        CTs=tris.CTs;
        nodeCTs=CT.nodes; areaChildCTs=CT.areaChild;iPs=CT.iPs;
        fprintf("%fs\n",toc-time);
    end
    
    
    %% Nodes
    msh.nNew = nodeCTs.n; % number of new nodes (nodes around the interface)
    msh.nStd = size(points,2); % number of standard nodes
    msh.ndof = msh.nNew + msh.nStd; % number of dofs
    msh.newNodes = getNewNodes(nodeCTs.all,msh.nStd); % vector contaning new numbering of nodes around interface, column
    msh.node = getNodes(tris,nodeCTs,msh,phi,pa); % get all nodes
    
    
    %% Boundary nodes
    [iN,bN] = getibNodes(msh);
    b3Nodes = bN.e3;
end