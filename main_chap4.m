%% =======================================================================
% This file is used to find a very simple SYSTEM of equations proposed by
%   Thi for the ARTICLE 1 & CHAPTER 4
% We need to solve w,v first and then u.
% ------------------------------------------------------------------------
% PURPOSE: Verify nxfem with nonlinear system.
% Related files: 
%   - ffpp-article1.edp: for the code in FreeFem++ (fitted mesh)
%   - ffpp-article1-cr.edp: for the convergence rate
%   - model_article1.m: contains all info for the model
%   - defG.m
% ------------------------------------------------------------------------
% UPDATE 10/11/18: come back after a very long time, edit from 2 old files
% main_article1.m and main_article1_each.m to only 1 file like main_nxfem.
% ------------------------------------------------------------------------


%% =======================================================================
% DOMAIN: [-1,1]X[-1,1], interface: circle ((0,0);r0) inside domain
%-------------------------------------------------------------------------
% MODELS:
% r=sqrt(x^2+y^2)-r0 the interface
%-------------------------------------------------------------------------
% Equation of w:
% -nabla(alpha*nabla(w)) = fw   in Omega
% [w]=[alpha*nabla_n(w)]=0     on Gamma
% w=wex   on \partial\Omega
%-------------------------------------------------------------------------
% Equation of v:
% -nabla(beta*nabla(v)) - lam*v*g(u) = fv in Omega
% v=nabla_n(v)=0 on Gamma  <== SPECIAL HERE!!!
% v=vex on \partial\Omega
%-------------------------------------------------------------------------
% Equation of u: u = w - beta/(alpha*lambda)*v
%-------------------------------------------------------------------------
% EXACT SOLUTION:
% uex = r^2/alp1 if r<=r0
%       (r^2-r0^2)/alp2 + r0^2/alp1 if r>r0
% vex = (r^2-r0^2)^2/bet1 if r<=r0  <== SPECIAL HERE!!!
%       0 if r>r0                   <== SPECIAL HERE!!!
% wex = uex + bet/(lam*alp)*vex
%-------------------------------------------------------------------------
% RHS:
% fu = -4+vex*g(uex) = -4 + (r^2-r0^2)^2/bet1*g(r^2/alp1) if r<=r0
%                      -4 if r>r0
% fv = 8(r0^2-2r^2)-lam*vex*g(uex) if r<=r0
%      0 if r<r0
% fw = fu + 1/lam*fv = -4 + 8/lam*(r0^2-2r^2) if r<=r0
%                      -4 if r>r0
%=========================================================================


%% add path of functions
addpath(genpath('func')); % add all necessary functions



%% Fixed parameters
pa.degP1D = 3; % Gaussian quadrature points in 1D (for polinomial functions)
pa.degP2D = 4; % Gaussian quadrature points in 2D (for polinomial functions)
pa.degN = 8; % Gaussian quadrature points in 2D (for non-polynomial functions)
% degree-#OfPoints : 1-1, 2-3, 3-4, 4-6, 5-7, 6-12, 7-13, 8-16, 9-19, 10-25, 11-27, 12-33
pa.tol = eps(1e3); % tolerance, 1e-14



%% Model
model = model_article1;
GeoDom = model.domain(); % domain



%% SETTINGS
findCR = 1; % wanna find the CR?
    %     numSegCR = [16, 32, 64, 128]; % only works with findCR=1
    numSegCR = [36, 56, 86, 126];
    showPlotCR = 1; % show plot of convergence (for findCR=1)
    
numSegPlot = 51; % only for plotting, findCR=0
savePlot = 0; % 1 = export figures to files (and without plotting them)
showPlot = 1; % wanna plot or not the solution? (JUST FOR savePlot=0)
    nf = 0; % counter of figures (plot each plot in a separated figure)
    
pa.smallCut = 0; % ignore small-support basis (1=ignore,0=no)
    pa.tH = 1e2; % to find the small support using (20) or (21) in arnold 2008
    
reguMesh = 1; % regular or irregular mesh?

useNewton = 0; % use Newton method for solving v?

% Penalty (goes with \int [][])
cpW.lamH = 1e11; % penalty coefficient for w
cpV.lamH = 1e8; % penalty coefficient for v

% ghost penalty
pa.useGP = 0; % wanna use ghost penalty term?
    pa.gam1 = 1e-6; % parameter for 1st term
    pa.gam2 = 1e-6 ; % parameter for 2nd term



%% Model's parameters
cpW.kk1 = 1; cpW.kk2 = 100;
cpV.kk1 = 0.5; cpV.kk2 = 100;
% (cpV.kk2 doesn't take affect because v=0 in Omg2)
pa.r0 = 0.6; % interface
pa.lamSys = 1; % coef lam in system settings



%% Dependent parameters
if findCR == 1
    showPlot = 0; % don't plot the solutions
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
    bNodes = bN.all; iNodes=iN.all;
    
    
    %% Exact solution in stdFEM
    % exSol_i = exSol(x_i)
    % w
    defExSol = model.defWex;
    wExStd = exInStd(defExSol,msh,pa);
    % u
    defExSol = model.defUex;
    uExStd = exInStd(defExSol,msh,pa);
    % v
    defExSol = model.defVex;
    vExStd = exInStd(defExSol,msh,pa);
    
    
    %% Exact solution in NXFEM
    % wExNX_i = wExStd_i for i is node of mesh
    % wExNX_k(i) = wExStd_i for i in msh.node.CT.all
    % w
    
    % u
    
    % v
    
    
end


%% =======================================================================
% SETTINGS
%=========================================================================
nStep = 4; % number of steps
% nStep = 1; % 1 step, for plotting
nSeg = 15;
useFilesNotPlot = 1; % used to run on server (without plotting pop-up)


if nStep>1
    numSeg = [21 41 61 91];
    % numSeg = [81 91 101 121];
else
    numSeg(1)=nSeg; 
end
for z = nStep:-1:1
    [hTmax(z),errW(z),errV(z)] = main_article1_each(numSeg(z));
end

%% finding order of CR
% -----------------------------------------------------------------------
if nStep>1 
    tx=zeros(1,nStep); 
    tyV=zeros(2,nStep);
    tyW=zeros(2,nStep);
    for i=1:nStep
       tx(i) = hTmax(i);
       
       tyV(1,i) = errV(i).L2;
       tyV(2,i) = errV(i).ENorm;
       
       tyW(1,i) = errW(i).L2;
       tyW(2,i) = errW(i).ENorm;
    end
    tmp = polyfit(log(tx),log(tyW(1,:)),1);
    orderW.L2 = tmp(1);
    tmp = polyfit(log(tx),log(tyW(2,:)),1);
    orderW.ENorm = tmp(1);

    tmp = polyfit(log(tx),log(tyV(1,:)),1);
    orderV.L2 = tmp(1);
    tmp = polyfit(log(tx),log(tyV(2,:)),1);
    orderV.ENorm = tmp(1);

    % Display results
    % --------------------------------------
    fprintf('L2 order of V: %0.15f\n',orderV.L2);
    fprintf('ENorm order of V: %0.15f\n',orderV.ENorm);
    
    fprintf('L2 order of W: %0.15f\n',orderW.L2);
    fprintf('ENorm order of W: %0.15f\n',orderW.ENorm);
end


%% Export to files/figures
% ------------------------------------------------------------------------
% See the file results\article1\...
% THIS IS NOT WORKING FOR A MOMENT!
if (useFilesNotPlot) && (nStep>1)
    fv=figure('visible','off');
    plot(log(tx),log(tyV(1,:)),'-r.',...
         log(tx),log(tyV(2,:)),'-b.');
    legend('L2V','ENormV');
    xlabel('log(h)'); 
    ylabel('log(errorV)');
    print -djpeg results\article1\errV.jpg;
    close(fv);
    
    fw=figure('visible','off');
    plot(log(tx),log(tyW(1,:)),'-r.',...
         log(tx),log(tyW(2,:)),'-b.');
    legend('L2W','ENormW');
    xlabel('log(h)'); 
    ylabel('log(errorW)');
    print -djpeg results\article1\errW.jpg;
    close(fw);
    
    % compute CR at each step
    crV = zeros(2,nStep);
    crW = zeros(2,nStep);
    for i=2:nStep
       crV(1,i) = log(tyV(1,i)/tyV(1,i-1))/log(tx(i)/tx(i-1)); % L2
       crV(2,i) = log(tyV(2,i)/tyV(2,i-1))/log(tx(i)/tx(i-1)); % ENorm
       
       crW(1,i) = log(tyW(1,i)/tyW(1,i-1))/log(tx(i)/tx(i-1)); % L2
       crW(2,i) = log(tyW(2,i)/tyW(2,i-1))/log(tx(i)/tx(i-1)); % ENorm
    end

    % see file results\article1\errorsV.txt
    errFileMatV = [tx;tyV(1,:);crV(1,:);tyV(2,:);crV(2,:)];
    fileID = fopen('results\article1\errorsV.txt','w');
    fprintf(fileID,'%7s %12s %6s %12s %6s\n','h','L2V','CR','ENormV','CR');
    fprintf(fileID,'%6.5f %12.8f %6.2f %12.8f %6.2f \n',errFileMatV);
    fclose(fileID);
    
    errFileMatW = [tx;tyW(1,:);crW(1,:);tyW(2,:);crW(2,:)];
    fileID = fopen('results\article1\errorsW.txt','w');
    fprintf(fileID,'%7s %12s %6s %12s %6s\n','h','L2w','CR','ENormW','CR');
    fprintf(fileID,'%6.5f %12.8f %6.2f %12.8f %6.2f \n',errFileMatW);
    fclose(fileID);
elseif nStep > 1
%% Plot the CR
    figure(1);
    plot(log(tx),log(tyV(1,:)),'-r.',...
         log(tx),log(tyV(2,:)),'-b.');
    legend('L2','ENorm');
    xlabel('log(h)'); 
    ylabel('log(error-V)');

    figure(2);
    plot(log(tx),log(tyW(1,:)),'-r.',...
         log(tx),log(tyW(2,:)),'-b.');
    legend('L2','ENorm');
    xlabel('log(h)'); 
    ylabel('log(error-W)');
end

