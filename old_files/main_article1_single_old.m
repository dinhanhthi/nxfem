%% =======================================================================
% This file is used to find a very simple SYSTEM of equations proposed by
% Thi for the ARTICLE 1
% We need to solve w,v first and then u.
%-------------------------------------------------------------------------
% PURPOSE: Verify nxfem with nonlinear system.
% Related files: 
%   - ffpp-article1.edp: for the code in FreeFem++ (fitted mesh)
%   - model_article1.m: contains all info for the model

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
model = model_article1; % choose model. cf. file model_sys_linda.m
% typeG: modify directly in defG.m

%-------------------------------------------------------------------------
% Settings
%-------------------------------------------------------------------------
wannaPlot = 1; % wanna plot?

%-------------------------------------------------------------------------
% Deleting small cut
%-------------------------------------------------------------------------
pa.smallCut = 1; % ignore small-support basis (1=ignore,0=no)
pa.tH = 10; % to find the small support using (20) or (21) in arnold 2008

%-------------------------------------------------------------------------
% Penalty parameters (\gam\int [w][varphi])
%-------------------------------------------------------------------------
pa.lamHw = 1e5; % penalty coefficient for w
pa.lamHv = 1e5; % penalty coefficient for v

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
nSeg = 45; % ONLY FOR nStep=1;

%-------------------------------------------------------------------------
% Model parameters
%-------------------------------------------------------------------------
pa.alp1 = 1; pa.alp2 = 100; % diff coeff for eqn w
pa.bet1 = 0.5; pa.bet2 = 100; % diff coeff for eqn v
% (bet2 doesn't take affect because v=0 in Omg2)
pa.r0 = 0.6; % interface
pa.lamSys = 1; % coef lam in system settings
pa.vgam1 = 1e5; % penalty term for equation of v (it's not ghost penalty)
pa.vgam2 = 1e3; % in Omg2



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
    [tri2del,t2Omg1,t2Omg2] = findSmallPhi(msh,CTs,CT,pa,phi);
    % If there are small-cut triangles, remove them from CTs!!
    if ~isempty(tri2del)
        nCTs = size(CTs,2); % number of OLD cut triangles
        % get NEW not-cut triangles
        if ~isempty(t2Omg1)
            NCTs1 = [NCTs1,CTs(:,t2Omg1)]; % add more triangles to NCTs1
            tris.NCTs1=NCTs1;
        end
        if ~isempty(t2Omg2)
            NCTs2 = [NCTs2,CTs(:,t2Omg2)]; % add more triangles to NCTs2
            tris.NCTs2=NCTs2;
        end
        % get NEW cut triangles
        CTs = CTs(:,setdiff(1:nCTs,tri2del));
        tris.CTs=CTs;
        % find again all information
        clear nodeCTs typeCTs areaChildCTs iPs uNVCTs areaCTs; % in case
        CT = getInfoCTs(CTs,phi,msh,pa);
        nodeCTs=CT.nodes; areaChildCTs=CT.areaChild;iPs=CT.iPs;
    end
end



%% =======================================================================
% NODES
%=========================================================================
msh.nNew = nodeCTs.n; % number of new nodes (nodes around the interface)
msh.nStd = size(points,2); % number of standard nodes
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

% w
%-------------------------------------------------------------------------
exSolNXw = interSTD2NX(exSolw,msh); % column array

% v
%-------------------------------------------------------------------------
exSolNXv = interSTD2NX(exSolv,msh); % column array



%% =======================================================================
% CONTROL PARAMETERS
% depend on mesh and different for tw and u
% in child-functions, it's the variable "cp"
%=========================================================================
kapW = model.kapW(areaChildCTs,pa);
cpW.kap1 = kapW.kap1; cpW.kap2 = kapW.kap2; % kappa_i
cpW.lambda = model.lamW(areaChildCTs,pa); % penalty coef (not ghost penalty)
cpW.kk1 = pa.alp1; cpW.kk2 = pa.alp2;    % diff coef

kapV = model.kapV(areaChildCTs,pa);
cpV.kap1 = kapV.kap1; cpV.kap2 = kapV.kap2; % kappa_i
cpV.lambda = model.lamV(areaChildCTs,pa); % penalty coef (not ghost penalty)
cpV.kk1 = pa.bet1; cpV.kk2 = pa.bet2;    % diff coef



%% =======================================================================
% SOLVING W
%=========================================================================

%-------------------------------------------------------------------------
% Stiffness matrix (all nodes including nodes on boundary)
%-------------------------------------------------------------------------
Aw = getGMGG(tris,phi,CT,msh,pa,cpW);

%-------------------------------------------------------------------------
% Load vector (all nodes including nodes on boundary)
%-------------------------------------------------------------------------
defFw = model.defFw;
Fw = getLf(tris,CT,msh,pa,defFw);

%-------------------------------------------------------------------------
% BCs
%-------------------------------------------------------------------------
numSolw = zeros(msh.ndof,1); % column-array
typeBC = model.bcW(); % get type of BCs
switch typeBC
    case 1 % u=o on whole boundary
        numSolw(bNodes) = 0;
    case 2 % u=uex on whole boundary
        numSolw(bNodes) = exSolNXw(bNodes);
end

%-------------------------------------------------------------------------
% Solving w
%-------------------------------------------------------------------------
Fw = Fw - Aw*numSolw; % modification of F
% LU factorization
numSolw(iNodes) = Aw(iNodes,iNodes)\Fw(iNodes); % don't care nodes on boundary
% numSoltw(iNodes) = gmres(Atw(iNodes,iNodes),Ftw(iNodes)); % GMRES factorization
sol.wh = numSolw;



%% =======================================================================
% SOLVING V
%=========================================================================
useNewton = 1; % wanna use Newton method or not?

% initial solution
% ------------------------------------------------------------------------
vold = zeros(msh.ndof,1);
% numSolu = ones(msh.ndof,1);
% numSolu = exSolNXu-0.05;

% tolerance
% ------------------------------------------------------------------------
itol = 1e-6;
delL2 = getErrL2(vold - exSolNXv,tris,CT,msh,pa);
Vip1L2 = getErrL2(exSolNXv,tris,CT,msh,pa);
difv = delL2/Vip1L2; % |del|_L2/|v_i+1|_L2fprintf('difv: %0.7f\n',difv);
fprintf('difv: %0.7f\n',difv);

imax = 50; % maximum number of steps
step = 0;
defFv = model.defFv;
typeBC = model.bcV(); % get type of BCs
wS = getUold(numSolw,msh); % w in each subdomain
while (difv > itol) && (step<imax)
    step = step+1;
    
    % analyze numSolu into each subdomain
    % --------------------------------------------------------------------
    voldEach = getUold(vold,msh);
    
    if ~useNewton % don't wanna use Newton method
        
        % global matrix
        % -----------------------------------------------------------------
        Av = getGMvAA(tris,phi,voldEach,wS,CT,msh,pa,cpV,cpW);
        
        % load vector
        % -----------------------------------------------------------------
        Fv = getLf(tris,CT,msh,pa,defFv);
        
        % bc
        % -----------------------------------------------------------------
        vnew = zeros(msh.ndof,1); % zero initial uh for each step
        switch typeBC
            case 1 % u=o on whole boundary
                vnew(bNodes) = 0;
            case 2 % u=uex on whole boundary
                vnew(bNodes) = exSolNXv(bNodes);
        end
        
        % solving for u
        % -----------------------------------------------------------------
        Fv = Fv - Av*vnew; % modification of F
        % LU factorization
        vnew(iNodes) = Av(iNodes,iNodes)\Fv(iNodes); % don't care nodes on boundary
        % uh(iNodes) = gmres(Au(iNodes,iNodes),Fu(iNodes)); % GMRES factorization
        
        del = vnew - vold;
        vold = vnew; % update for the next step
        
    else % use Newton method (DF(u)*del = F(u), solve for del)
        
        % DF(u)*del
        % -----------------------------------------------------------------
        Adel = getGMvAANewton(tris,phi,voldEach,wS,CT,msh,pa,cpV,cpW);
               
        % F(u)
        % -----------------------------------------------------------------
        
        Av = getGMvAA(tris,phi,voldEach,wS,CT,msh,pa,cpV,cpW); 
                    % like Av in normal iterative method
        Fv = getLf(tris,CT,msh,pa,defFv);
                    % like Fu in normal iterative method
        Fdel = Av*vold - Fv;
        
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
        
        vold(bNodes) = exSolNXv(bNodes); % u=uex on bc
        vold(iNodes) = vold(iNodes) - del(iNodes); % u_i+1 = u_i - del, update for the next step
    end
    
    delL2 = getErrL2(del,tris,CT,msh,pa);
    Vip1L2 = getErrL2(vold,tris,CT,msh,pa);
%     difu = delL2/Uip1L2; % |del|_L2/|u_i+1|_L2
    difv = delL2;
    fprintf('difv: %0.18f\n',difv);
end
sol.vh = vold;



%% =======================================================================
% SOLVING U
% u = w - bet/(alp*lamSys)*v
%=========================================================================




%% ========================================================
% NUMERICAL SOLUTION 
% in STANDARD FEM (just for plotting)
% VhSol_i = numSol_i for i in nodesOmg1NotOnGam or nodesOmg2NotAroundGam
% VhSol_i = numSol_k(i) for i in nodesOmg2CT
% VhSol_i = numSol_i+numSol_k(i) for i in nodesOnGam
% =========================================================

% w
%-------------------------------------------------------------------------
VhSolw = interNX2STD(numSolw,msh);
sol.wVh = VhSolw; % add to var

% v
%-------------------------------------------------------------------------
VhSolv = interNX2STD(vold,msh);
sol.vVh = VhSolv; % add to var



%% =======================================================================
% PLOTTING
%=========================================================================
if wannaPlot
    nf = 0;
    % nf = plotNXFEM(msh,iPs,nf,'eleLabel','off','nodeLabel','on'); % only mesh
    nf = plotNXFEM(msh,iPs,nf,sol.wVh,'withMesh',false,...
                    'title','wh','dim',3); % wh
    nf = plotNXFEM(msh,iPs,nf,sol.wex,'withMesh',false,'title','wex',...
                'dim',3,'export',false); % wex
    nf = plotNXFEM(msh,iPs,nf,sol.vVh,'withMesh',false,...
                    'title','vh','dim',3); % vh
    nf = plotNXFEM(msh,iPs,nf,sol.vex,'withMesh',false,'title','vex',...
                'dim',3,'export',false); % vex
end

