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
    itol = 1e-6;
    imax = 50; % number of iterative

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
    disp('Get mesh info...');
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
    
    
    disp('Exact solutions...');
    %% Exact solution in stdFEM
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
    wExNX = interSTD2NX(wExStd,msh); % column array
    % u
    uExNX = interSTD2NX(uExStd,msh); % column array
    % v
    vExNX = interSTD2NX(vExStd,msh); % column array
    
    
    %% Control paramaters
    hTCTs = msh.hT(CTs(5,:));
    % w
    kapW = model.kapW(cpW,CT,pa);
    cpW.kap1 = kapW.kap1; cpW.kap2 = kapW.kap2; % kappa_i
    cpW.lambda = model.lamW(cpW,hTCTs,CT,pa); % penalty coef (not ghost penalty)
    % v
    kapV = model.kapV(cpV,CT,pa);
    cpV.kap1 = kapV.kap1; cpV.kap2 = kapV.kap2; % kappa_i
    cpV.lambda = model.lamV(cpV,hTCTs,CT,pa); % penalty coef (not ghost penalty)
    
    
    %% SOLVING W
    fprintf('Solving w... ');tic;time=0;
    Aw = getGMGG(tris,phi,CT,msh,pa,cpW);
    
    defFw = model.defFw;
    Fw = getLf(tris,CT,msh,pa,defFw);
    
    whNX = zeros(msh.ndof,1); % column-array
    typeBC = model.bcW(); % get type of BCs
    switch typeBC
        case 1 % u=o on whole boundary
            whNX(bNodes) = 0;
        case 2 % u=uex on whole boundary
            whNX(bNodes) = wExNX(bNodes);
    end
    
    Fw = Fw - Aw*whNX;
    whNX(iNodes) = Aw(iNodes,iNodes)\Fw(iNodes); % don't care nodes on boundary
    % whNX(iNodes) = gmres(Aw(iNodes,iNodes),Fw(iNodes)); % GMRES factorization
    
    fprintf("%fs\n",toc-time);
    
    
    %% SOLVING V
    fprintf('Solving v... ');tic;time=0;
    
    % initial v
    difu = 100; % initial, harmless
    step = 0;
    defFv = model.defFv;
    typeBC = model.bcV(); % get type of BCs
    defG = defGu;
    wS = getUold(whNX,msh); % w in each subdomain
    
    if ~useNewton
        fprintf('Normal fixed point method)\n');
    else
        fprintf('Using Newton method)\n');
    end
    
    while (difu > itol) && (step<imax)
        step = step+1;
        voldEach = getUold(vold,msh);
        if ~useNewton % don't use Newton
            Av = getGMvAA(tris,phi,voldEach,wS,CT,msh,pa,cpV,cpW);
            Fv = getLf(tris,CT,msh,pa,defFv);
            
            vnew = zeros(msh.ndof,1); % zero initial uh for each step
            switch typeBC
                case 1 % u=o on whole boundary
                    vnew(bNodes) = 0;
                case 2 % u=uex on whole boundary
                    vnew(bNodes) = vExNX(bNodes);
            end
            
            Fv = Fv - Av*vnew; % modification of F
            vnew(iNodes) = Av(iNodes,iNodes)\Fv(iNodes); % don't care nodes on boundary
            % vnew(iNodes) = gmres(Av(iNodes,iNodes),Fv(iNodes)); % GMRES factorization
            
            del = vnew - vold;
            vold = vnew;
        else
            % DF(u)
            Adel = getGMvAANewton(tris,phi,voldEach,wS,CT,msh,pa,cpV,cpW);
            % F(u)
            Av = getGMvAA(tris,phi,voldEach,wS,CT,msh,pa,cpV,cpW); % like Av in normal iterative method
            Fv = getLf(tris,CT,msh,pa,defFv); % like Fu in normal iterative method
            Fdel = Av*vold - Fv;
            
            del = zeros(msh.ndof,1); % column-array
            del(bNodes) = 0; % always
            
            Fdel = Fdel - Adel*del;
            del(iNodes) = Adel(iNodes,iNodes)\Fdel(iNodes); % don't care nodes on boundary
            % del(iNodes) = gmres(Adel(iNodes,iNodes),Fdel(iNodes)); % GMRES factorization
            
            vold(bNodes) = vExNX(bNodes); % v=vex on bc
            vold(iNodes) = vold(iNodes) - del(iNodes);
        end
        
        delL2 = getNormL2fhNX(del,tris,CT,msh,pa);
        Uip1L2 = getNormL2fhNX(vold,tris,CT,msh,pa);
        difu = delL2/Uip1L2; % |del|_L2/|u_i+1|_L2
        fprintf('___difu: %0.18f\n',difu);
    end
    
    vhNX = vold;
    
    
    
    %% GET U
    uhNX = getUh(whNX,vhNX,pa.bet1/(pa.lamSys*pa.alp1),...
                    pa.bet2/(pa.lamSys*pa.alp2),msh);
                
    
                
    %% ERRORS
    
end



