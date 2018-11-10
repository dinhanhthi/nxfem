%% ========================================================
% INTRODUCTION
% =========================================================
% NXFEM based on Hansbo's paper.
% This file is used for many models have a very simple form
%   -grad(kgrad u) = f
%   with a given exact solution
% =========================================================

%% add path of functions
addpath(genpath('func')); % add all necessary functions



%% Fixed parameters
pa.degP1D = 3; % Gaussian quadrature points in 1D (for polinomial functions)
pa.degP2D = 4; % Gaussian quadrature points in 2D (for polinomial functions)
pa.degN = 8; % Gaussian quadrature points in 2D (for non-polynomial functions)
% degree-#OfPoints : 1-1, 2-3, 3-4, 4-6, 5-7, 6-12, 7-13, 8-16, 9-19, 10-25, 11-27, 12-33
pa.tol = eps(1e3); % tolerance, 1e-14



%% model 
mdl = 2; 
switch mdl
    case 1 % sinha (cf. model_sinha.m)
        model=model_sinha; % file sinha.m
        cp.kk1 = 1; cp.kk2 = 0.5; % diffusion coefficients (kk1 must be = 2kk2)
        pa.xi = 1.0; % DON'T CHANGE!!!
    case 2 % barrau_x (p.44) STRAIGHT LINE [-1,1]x[-1,1]
        model=model_barrau_x; % file barrau_x.m
        cp.kk1 = 1; cp.kk2 = 1000; % for barrau case
        pa.xi = 0.31; % only used for barrau's model, between (-1,1)
    case 3 % barrau_r (p.47) PART OF CIRCLE [0,1]x[0,1]
        model=model_barrau_r; % file barrau_r.m
        cp.kk1 = 1; cp.kk2 = 1000; % for barrau case
        pa.xi = 0.75; % only used for barrau's model, between (0,1)
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
pa.kk1 = cp.kk1; pa.kk2 = cp.kk2; % just for defF



%% SETTINGS
findCR = 1; % wanna find convergence rate? 1 or 0
%     numSegCR = [16, 32, 64, 128]; % only works with findCR=1
    numSegCR = [36, 56, 86, 126];
    showPlotCR = 1; % show plot of convergence (for findCR=1)
numSegPlot = 51; % only for plotting, findCR=0
savePlot = 0; % 1 = export figures to files (and without plotting them)
    showPlot = 1; % wanna plot or not the solution? (JUST FOR savePlot=0)
    nf = 0; % counter of figures (plot each plot in a separated figure)
pa.smallCut = 0; % ignore small-support basis (1=ignore,0=no)
    pa.tH = 100; % to find the small support using (20) or (21) in arnold 2008
reguMesh = 1; % regular or irregular mesh?
    
cp.lamH = 1e4; % penalty coefficient

% ghost penalty
pa.useGP = 1; % wanna use ghost penalty term?
    pa.gam1 = 1e-9; % parameter for 1st term
    pa.gam2 = 1e-9; % parameter for 2nd term

    
    
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
        clear CTs NCTs NCTs2 nodeCTs areaChildCTs iPs; % just in case
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
    bNodes = bN.all; iNodes=iN.all;
    
    
    %% Exact solution in stdFEM
    % exSol_i = exSol(x_i)
    defExSol = model.defExSol; % function handle
    uExStd = exInStd(defExSol,msh,pa);


    %% Exact solution in NXFEM
    % wExNX_i = wExStd_i for i is node of mesh
    % wExNX_k(i) = wExStd_i for i in msh.node.CT.all
    uExNX = interSTD2NX(uExStd,msh); % column array


    %% Control paramaters
    % depend on mesh
    hTCTs = msh.hT(CTs(5,:));
    kap = model.kap(cp,CT,pa);
    cp.kap1 = kap.kap1; cp.kap2 = kap.kap2; % kappa_i
    cp.lambda = model.lam(cp,hTCTs,CT,pa); % penalty coef (not ghost penalty)

    
    fprintf('Solving u... ');tic;time=0;
    %% Stiffness
    % on all nodes including nodes on the boundary
    A = getGMGG(tris,phi,CT,msh,pa,cp);


    %% Load vector
    defF = model.defF;
    F = getLf(msh,pa,tris,CT,defF);


    %% Applying boundary conditions
    uhNX = zeros(msh.ndof,1); % column-array
    typeBC = model.bc(); % get type of BCs
    switch typeBC
        case 1 % u=o on whole boundary
            uhNX(bNodes) = 0;
        case 2 % u=uex on whole boundary
            uhNX(bNodes) = uExNX(bNodes);
    end


    %% Get solution
    F = F - A*uhNX;
    % LU factorization
    uhNX(iNodes) = A(iNodes,iNodes)\F(iNodes); % don't care nodes on boundary
    % uhNX(iNodes) = gmres(A(iNodes,iNodes),F(iNodes)); % GMRES factorization
    fprintf("%fs\n",toc-time);
    

    %% num sol in stdFEM
    % VhSol_i = numSol_i for i in nodesOmg1NotOnGam or nodesOmg2NotAroundGam
    % VhSol_i = numSol_k(i) for i in nodesOmg2CT
    % VhSol_i = numSol_i+numSol_k(i) for i in msh.node.onG
    uhStd = interNX2STD(uhNX,msh);
    
    
    %% Get errors
    disp('Get errors...');
    eU = uhNX - uExNX; % wex-wh
    
    fprintf('__L2... ');tic;time=0;
    L2 = getNormL2uhuexNX(msh,pa,tris,CT,uhNX,defExSol,defPhi);
    fprintf("%fs\n",toc-time);
    
    % L2h: testing
    fprintf('__L2h... ');tic;time=0;
    L2h = getNormL2fhNX(eU,tris,CT,msh,pa);
    fprintf("%fs\n",toc-time);
    
    % IMPORTANT (LATER): Need to write functions compute uex (function
    %   handle) and uh, not eU as below!!!!!
%     fprintf('__L2Gh... ');tic;time=0;
%     L2G = getNormL2GfhNX(eU,tris,areaChildCTs,msh,cp); % ||kgrad||_L2
%     jumU = getNormJump1p2(eU,CTs,iPs,msh,pa);
%     avegnU = getNormAveGn(eU,CTs,iPs,msh,cp); % ||{kgran w}||_{-1/2}
%     ENorm = L2G^2 + jumU^2+ avegnU^2;
%     ENorm = sqrt(ENorm);
%     fprintf("%fs\n",toc-time);
    
    
    %% For finding CR
    hTarray(z) = msh.hTmax;
    nDOFsArray(z) = msh.ndof; % including new nodes
    errArray(1,z) = L2;
%     errArray(2,z) = ENorm;
    errArray(2,z) = L2h;
end



%% Display CR or plotting
if nStep>1
    
    %% find CR
    disp('Find the CR...');
    tmp = polyfit(log(hTarray),log(errArray(1,:)),1);
    order.L2 = tmp(1);
    tmp = polyfit(log(hTarray),log(errArray(2,:)),1);
    order.ENorm = tmp(1);
    
    
    %% Save figure of CR
    if savePlot
        disp('Saving plot of CR...');
        
        %% Save figure of CR
        f = figure('visible','off');
        plot(log(hTarray),log(errArray(1,:)),'-r.',...
            log(hTarray),log(errArray(2,:)),'-b.');
        legend('L2','ENorm');
        xlabel('log(h)'); 
        ylabel('log(error)');
        print -djpeg results\main_nxfem\err_L2_ENorm.jpg;
        close(f);
        
        %% save errors as file
        % see file results\main_nxfem\errors.txt
        cr = zeros(2,nStep);
        for i=2:nStep
           cr(1,i) = log(errArray(1,i)/errArray(1,i-1))/log(hTarray(i)/hTarray(i-1)); % L2
           cr(2,i) = log(errArray(2,i)/errArray(2,i-1))/log(hTarray(i)/hTarray(i-1)); % ENorm
        end
        errFileMat = [hTarray;errArray(1,:);cr(1,:);errArray(2,:);cr(2,:)];
        fileID = fopen('results\main_nxfem\err_L2_ENorm.txt','w');
        fprintf(fileID,'Model: %2.0f \n',mdl);
        fprintf(fileID,'%7s %12s %6s %12s %6s\n','h','L2','CR','ENorm','CR');
        fprintf(fileID,'%6.5f & %12.8f & %6.2f & %12.8f & %6.2f \n',errFileMat);
        fclose(fileID);
    end
    
    
    %% Show plots CR
    if showPlotCR
        disp('Showing plot of CR...');
        plot(log(hTarray),log(errArray(1,:)),'-r.',...
            log(hTarray),log(errArray(2,:)),'-b.');
        legend('L2','ENorm');
        xlabel('log(h)'); 
        ylabel('log(error)');
    end
    
    %% Display CR
    fprintf('regular: %d,\t', reguMesh);
    fprintf('max numSeg: %0.0f,\t', numSeg(nStep));
    fprintf('max DOFs: %0.0f\n', nDOFsArray(nStep));
    fprintf('L2 order: %0.15f\n',order.L2);
%     fprintf('ENorm order: %0.15f\n',order.ENorm);
    fprintf('L2h order: %0.15f\n',order.ENorm);
    
else % Just for plotting, don't wanna fine CR
    
    fprintf('numSeg: %0.0f,\t',numSeg(nStep));
    fprintf('DOFs: %0.0f\n',nDOFsArray(nStep));
    fprintf('L2 err: %0.18f\n',L2);
    fprintf('ENorm err: %0.18f\n',ENorm);
    
    if showPlot
%         nf = plotNXFEM(msh,pa,phi,iPs,nf,'eleLabel','off','nodeLabel','off'); % only mesh
        nf = plotNXFEM(msh,pa,phi,iPs,nf,uhStd,'withMesh',false,'title','uh','dim',3); % uh
        nf = plotNXFEM(msh,pa,phi,iPs,nf,uExStd,'withMesh',false,'title','uex','dim',3); % uex
    end
end
