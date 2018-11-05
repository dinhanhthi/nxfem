%% =======================================================================
% This file is used to find a very simple SYSTEM of equations proposed by
% Linda. Note that this system is not like the system given in article 1
% because we don't have v=nabla_n v=0 on Gamma and we will solve w,u first
% and then v while in article 1, we need to solve w,v first and then u.
%-------------------------------------------------------------------------
% PURPOSE: Verify nxfem with nonlinear system.
% cf. linda-simple-system-ffpp.edp for the code in FreeFem++ (fitted mesh)
%-------------------------------------------------------------------------
% OLD FILE: from modifying getGMgPP for a more general case used in chopp06combine
%-------------------------------------------------------------------------
% DOMAIN: [0,1]X[0,1], interface: circle ((0.5,0.5);r0) inside domain
%-------------------------------------------------------------------------
% MODELS:
% Note that we find tw instead of w (w=beta*tw)
% r=sqrt((x-.5)^2+(y-.5)^2)
% Equation of tw:
% -nabla(beta*nabla(tw)) = -16r^2   in Omega
% [tw]=[beta*nabla_n(tw)]=0     on Gamma
% tw=twex   on \partial\Omega
% Equation of u:
% -nabla(alpha*nabla(u)) + beta*v*g(u) = -4+(wex-alpha*uex)*g(uex) in Omega
% [u]=[alpha*nabla_n(u)]=0 on Gamma
% u=uex on \partial\Omega
% Equation of v: v = 1/beta*(w-alpha*u)
%-------------------------------------------------------------------------
% EXACT SOLUTION:
% twex = r^4/beta1  if r<=r0
%        r^4/beta2 - r0^4/beta2 + r0^4/beta1  if r>r0
% uex = r^2/alpha1  if r<=r0
%       r^2/alpha2 - r0^2/alpha2 + r0^2/alpha1  if r>r0
%=========================================================================

addpath(genpath('func')); % add all necessary functions

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
% degree-#OfPoints : 1-1, 2-3, 3-4, 4-6, 5-7, 6-12,
%                    7-13, 8-16, 9-19, 10-25, 11-27, 12-33
pa.tol = eps(1e3); % tolerance, 1e-14
model = model_sys_linda; % choose model. cf. file model_sys_linda.m
% typeG: modify directly in defG.m

%-------------------------------------------------------------------------
% Settings
%-------------------------------------------------------------------------
wannaPlot = 0; % wanna plot?

%-------------------------------------------------------------------------
% Deleting small cut
%-------------------------------------------------------------------------
pa.smallCut = 1; % ignore small-support basis (1=ignore,0=no)
pa.tH = 10; % to find the small support using (20) or (21) in arnold 2008

%-------------------------------------------------------------------------
% Penalty parameters (\gam\int [u][v])
%-------------------------------------------------------------------------
pa.lamHw = 1e5; % penalty coefficient for tw
pa.lamHu = 1e5; % penalty coefficient for u

%-------------------------------------------------------------------------
% Ghost penalty
%-------------------------------------------------------------------------
pa.useGP = 1; % wanna use ghost penalty term?
pa.gam1 = 1e-6; % parameter for 1st term
pa.gam2 = 1e-6 ; % parameter for 2nd term

%-------------------------------------------------------------------------
% Mesh settings
%-------------------------------------------------------------------------
pa.reguMesh = 0; % use regular mesh or not?
nSeg = 17; % ONLY FOR nStep=1;

%-------------------------------------------------------------------------
% Model parameters
%-------------------------------------------------------------------------
pa.alp1 = 1; pa.alp2 = 100; % diff coeff for eqn u
pa.bet1 = 1; pa.bet2 = 100; % diff coeff for eqn tw
pa.r0 = 0.3; % interface



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
msh.nStd = size(points,2); % number of standard nodes

%-------------------------------------------------------------------------
% Level set function
%-------------------------------------------------------------------------
phi = model.defPhi(x,y,pa); % 1 x number of points (row array)
phi(abs(phi)<pa.tol)=0; % find phi which are very small (~0) and set to 0



%% =======================================================================
% GET INFORMATION OF TRIANGLES
% The same for both equations of tw and u
%=========================================================================

%-------------------------------------------------------------------------
% Triangles
%-------------------------------------------------------------------------
tris = getTriangles(phi,msh,pa); % tris has 3 factors (structure var)
CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;

%-------------------------------------------------------------------------
% On cut triangles
%-------------------------------------------------------------------------
CT = getInfoCTs(CTs,phi,msh,pa); % CT has many factors (structure var)
nodeCTs=CT.nodes; areaChildCTs=CT.areaChild;iPs=CT.iPs;
        
%-------------------------------------------------------------------------
% Find small-cut triangles (idx in the OLD CTs)
%-------------------------------------------------------------------------
if pa.smallCut
    [tris,CT] = findSmallPhi_after(msh,pa,phi,tris,CT);
    clear CTs NCTs NCTs2 nodeCTs areaChildCTs iPs; % just in case
    CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
    nodeCTs=CT.nodes; areaChildCTs=CT.areaChild; iPs=CT.iPs;
end



%% =======================================================================
% NODES
%=========================================================================
msh.nNew = nodeCTs.n; % number of new nodes (nodes around the interface)
% msh.nStd = size(points,2); % number of standard nodes (already declared)
msh.ndof = msh.nNew + msh.nStd; % number of dofs
msh.newNodes = getNewNodes(nodeCTs.all,msh.nStd);
    % vector contaning new numbering of nodes around interface, column
msh.node = getNodes(tris,nodeCTs,msh,phi,pa); % get all nodes

%-------------------------------------------------------------------------
% boundary nodes and inner nodes
%-------------------------------------------------------------------------
[iN,bN] = getibNodes(msh);
bNodes = bN.all; iNodes=iN.all;



%% ========================================================
% EXACT SOLUTION
% in standard FEM space (just for plotting)
% exSol_i = exSol(x_i)
% =========================================================

% tw
%-------------------------------------------------------------------------
defExSol = model.defTWex;
exSoltw = exInStd(defExSol,msh,pa);
sol.twex = exSoltw; % add to var

% w
%-------------------------------------------------------------------------
defExSol = model.defWex;
exSolw = exInStd(defExSol,msh,pa);
sol.wex = exSolw; % add to var

% u
%-------------------------------------------------------------------------
defExSol = model.defUex;
exSolu = exInStd(defExSol,msh,pa);
sol.uex = exSolu; % add to var

% v
%-------------------------------------------------------------------------
defExSol = model.defVex;
exSolv = exInStd(defExSol,msh,pa);
sol.vex = exSolv; % add to var



%% ========================================================
% EXACT SOLUTION in NXFEM
% exSolNX_i = exSol_i for i is node of mesh
% exSolNX_k(i) = exSol_i for i in nodesAroundGam
% =========================================================

% tw
%-------------------------------------------------------------------------
exSolNXtw = interSTD2NX(exSoltw,msh); % column array

% u
%-------------------------------------------------------------------------
exSolNXu = interSTD2NX(exSolu,msh); % column array



%% =======================================================================
% CONTROL PARAMETERS
% depend on mesh and different for tw and u
% in child-functions, it's the variable "cp"
%=========================================================================
kapTW = model.kapTW(areaChildCTs,pa);
cpTW.kap1 = kapTW.kap1; cpTW.kap2 = kapTW.kap2; % kappa_i
cpTW.lambda = model.lamTW(areaChildCTs,pa); % penalty coef (not ghost penalty)
cpTW.kk1 = pa.bet1; cpTW.kk2 = pa.bet2;    % diff coef

kapU = model.kapU(areaChildCTs,pa);
cpU.kap1 = kapU.kap1; cpU.kap2 = kapU.kap2; % kappa_i
cpU.lambda = model.lamU(areaChildCTs,pa); % penalty coef (not ghost penalty)
cpU.kk1 = pa.alp1; cpU.kk2 = pa.alp2;    % diff coef



%% =======================================================================
% SOLVING TW
% don't forget tw=1/beta*w
%=========================================================================

%-------------------------------------------------------------------------
% Stiffness matrix (all nodes including nodes on boundary)
%-------------------------------------------------------------------------
Atw = getGMGG(tris,phi,CT,msh,pa,cpTW);

%-------------------------------------------------------------------------
% Load vector (all nodes including nodes on boundary)
%-------------------------------------------------------------------------
defFtw = model.defFtw;
Ftw = getLf(msh,pa,tris,CT,defFtw);



%-------------------------------------------------------------------------
% BCs
%-------------------------------------------------------------------------
numSoltw = zeros(msh.ndof,1); % column-array
typeBC = model.bcTW(); % get type of BCs
switch typeBC
    case 1 % u=o on whole boundary
        numSoltw(bNodes) = 0;
    case 2 % u=uex on whole boundary
        numSoltw(bNodes) = exSolNXtw(bNodes);
end

%-------------------------------------------------------------------------
% Solving tw
%-------------------------------------------------------------------------
Ftw = Ftw - Atw*numSoltw; % modification of F
% LU factorization
numSoltw(iNodes) = Atw(iNodes,iNodes)\Ftw(iNodes); % don't care nodes on boundary
% numSoltw(iNodes) = gmres(Atw(iNodes,iNodes),Ftw(iNodes)); % GMRES factorization
sol.twh = numSoltw;



%% =======================================================================
% W
% w=beta*tw
%=========================================================================
wS = getWsep(numSoltw,msh,pa.bet1,pa.bet2);



%% =======================================================================
% SOLVING U
% don't forget tw=1/beta*w
%=========================================================================
useNewton = 1; % wanna use Newton method or not?

% initial solution
% ------------------------------------------------------------------------
uold = zeros(msh.ndof,1);
% numSolu = ones(msh.ndof,1);
% numSolu = exSolNXu-0.05;

% tolerance
% ------------------------------------------------------------------------
itol = 1e-8;
delL2 = getNormL2(uold - exSolNXu,tris,CT,msh,pa);
Uip1L2 = getNormL2(exSolNXu,tris,CT,msh,pa);
difu = delL2/Uip1L2; % |del|_L2/|u_i+1|_L2fprintf('difu: %0.7f\n',difu);
fprintf('difu: %0.7f\n',difu);

imax = 50; % maximum number of steps
step = 0;
defG = defGu; % cf. defGu.m
defFu = @(xx,yy,pa,sub) model.defFu(xx,yy,pa,sub,defG.change);
typeBC = model.bcU(); % get type of BCs
while (difu > itol) && (step<imax)
    step = step+1;
    
    % analyze numSolu into each subdomain
    % --------------------------------------------------------------------
    uoldEach = getWsep(uold,msh,1,1);
    
    if ~useNewton % don't wanna use Newton method
        
        % global matrix
        % ----------------------------------------------------------------- 
        Au = getGMgPP(msh,pa,cpU,tris,CT,phi,uoldEach,defG);
        
        
        % load vector
        % -----------------------------------------------------------------          
        Fu = getLfwg(msh,pa,tris,CT,uoldEach,wS,defFu,defG);
        
        % bc
        % -----------------------------------------------------------------
        unew = zeros(msh.ndof,1); % zero initial uh for each step
        switch typeBC
            case 1 % u=o on whole boundary
                unew(bNodes) = 0;
            case 2 % u=uex on whole boundary
                unew(bNodes) = exSolNXu(bNodes);
        end
        
        % solving for u
        % -----------------------------------------------------------------
        Fu = Fu - Au*unew; % modification of F
        % LU factorization
        unew(iNodes) = Au(iNodes,iNodes)\Fu(iNodes); % don't care nodes on boundary
        % unew(iNodes) = gmres(Au(iNodes,iNodes),Fu(iNodes)); % GMRES factorization
        
        del = unew - uold;
        uold = unew; % update for the next step
        
    else % use Newton method (DF(u)*del = F(u), solve for del)
        
        % DF(u)*del
        % -----------------------------------------------------------------       
        Adel = getGMuNewton(msh,pa,cpU,tris,CT,phi,uoldEach,wS,defG);
                
        % F(u)
        % -----------------------------------------------------------------
        Au = getGMgPP(msh,pa,cpU,tris,CT,phi,uoldEach,defG);
                    % like Au in normal iterative method
        Fu = getLfwg(msh,pa,tris,CT,uoldEach,wS,defFu,defG);
                    % like Fu in normal iterative method
        Fdel = Au*uold - Fu;
        
        % bc for del
        % -----------------------------------------------------------------
        del = zeros(msh.ndof,1); % column-array
        del(bNodes) = 0; % always
        
        % solve for del
        % -----------------------------------------------------------------
        Fdel = Fdel - Adel*del; 
        % LU factorization
        del(iNodes) = Adel(iNodes,iNodes)\Fdel(iNodes); % don't care nodes on boundary
        % del(iNodes) = gmres(Adel(iNodes,iNodes),Fdek(iNodes)); % GMRES factorization
        
        uold(bNodes) = exSolNXu(bNodes); % u=uex on bc
        uold(iNodes) = uold(iNodes) - del(iNodes); % u_i+1 = u_i - del, update for the next step
    end
    
    delL2 = getNormL2(del,tris,CT,msh,pa);
    Uip1L2 = getNormL2(uold,tris,CT,msh,pa);
%     difu = delL2/Uip1L2; % |del|_L2/|u_i+1|_L2
    difu = delL2;
    fprintf('difu: %0.18f\n',difu);
end
sol.uh = uold;



%% ========================================================
% NUMERICAL SOLUTION 
% in STANDARD FEM (just for plotting)
% VhSol_i = numSol_i for i in nodesOmg1NotOnGam or nodesOmg2NotAroundGam
% VhSol_i = numSol_k(i) for i in nodesOmg2CT
% VhSol_i = numSol_i+numSol_k(i) for i in nodesOnGam
% =========================================================

% tw
%-------------------------------------------------------------------------
VhSoltw = interNX2STD(numSoltw,msh);
sol.twVh = VhSoltw; % add to var

% u
%-------------------------------------------------------------------------
VhSolu = interNX2STD(uold,msh);
sol.uVh = VhSolu; % add to var



%% =======================================================================
% PLOTTING
%=========================================================================
nf = 0;
% nf = plotNXFEM(msh,iPs,nf,'eleLabel','off','nodeLabel','on'); % only mesh
% nf = plotNXFEM(msh,iPs,nf,sol.twVh,'withMesh',false,...
%                 'title','twh','dim',3); % twh
% nf = plotNXFEM(msh,iPs,nf,sol.twex,'withMesh',false,'title','twex',...
%             'dim',3,'export',false); % twex
nf = plotNXFEM(msh,iPs,nf,sol.uVh,'withMesh',false,...
                'title','uh','dim',3); % uh
nf = plotNXFEM(msh,iPs,nf,sol.uex,'withMesh',false,'title','uex',...
            'dim',3,'export',false); % uex


