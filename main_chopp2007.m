%% =======================================================================
% This file is used to find a very simple BIOFILM proposed by Chopp 2007
% This is a SYSTEM of LINEAR equations u,v (u=substrate, v=potential)
% ------------------------------------------------------------------------
% REWRITE NOT FINISHED YET!
% ------------------------------------------------------------------------
% PURPOSE: Coding level set
% Related files: 
%   - ChoppSimpleKap.edp: for the code in FreeFem++
%   - model_chopp2007.m: contains all info for the model
% ------------------------------------------------------------------------


%% =======================================================================
% DOMAIN: [0,1]X[0,1], inital interface: half of circle ((1/2,0);r0)
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
model = model_chopp2007; % choose model. cf. file model_article1.m

%-------------------------------------------------------------------------
% Deleting small cut
%-------------------------------------------------------------------------
pa.smallCut = 1; % ignore small-support basis (1=ignore,0=no)
pa.tH = 10; % to find the small support using (20) or (21) in arnold 2008

%-------------------------------------------------------------------------
% Penalty parameters (\gam\int [w][varphi])
%-------------------------------------------------------------------------
pa.lamHu = 1e7; % penalty coefficient for u (substrate)
pa.lamHv = 1e7; % penalty coefficient for v (potential)

%-------------------------------------------------------------------------
% Ghost penalty
%-------------------------------------------------------------------------
pa.useGP = 1; % wanna use ghost penalty term?
pa.gam1 = 1e-7; % parameter for 1st term
pa.gam2 = 1e-7 ; % parameter for 2nd term

%-------------------------------------------------------------------------
% Mesh settings
%-------------------------------------------------------------------------
pa.reguMesh = 0; % use regular mesh or not?

%-------------------------------------------------------------------------
% Model parameters
%-------------------------------------------------------------------------
pa.alp1 = 120; pa.alp2 = 150; % diff coeff for eqn u
pa.bet1 = 1e6; pa.bet2 = 0; % this is not diff coef!
pa.r0 = 0.3; % interface
pa.mu1 = 3.6e6;
pa.mu2 = 0;
pa.bcu3 = 1e-3; % boundary condition for u on \pt\Omg_3



%% =======================================================================
% DOMAIN
%=========================================================================
GeoDom = model.domain(); % domain

%-------------------------------------------------------------------------
% Mesh setting up
%-------------------------------------------------------------------------
nSeg = 11; % mesh settings
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


%% =======================================================================
% EACH TIME STEP
%=========================================================================
pA = model.pa;
pA = pA();
t=0;
dt = pA.dt;
Tmax = pA.Tmax;
for ns = 1:maxStep
    
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
    b3Nodes = bN.e3; % node on \pt\Omg_3 (top)
    
    
    %% =======================================================================
    % CONTROL PARAMETERS
    % depend on mesh and different for w and u
    % in child-functions, it's the variable "cp"
    %=========================================================================
    kapU = model.kapU(areaChildCTs,pa);
    cpU.kap1 = kapU.kap1; cpU.kap2 = kapU.kap2; % kappa_i
    cpU.kk1 = pa.alp1; cpU.kk2 = pa.alp2;    % diff coef
    cpU.lambda = model.lamU(cpU,msh.hT,CTs,pa); % penalty coef (not ghost penalty)

    kapV = model.kapV(areaChildCTs,pa);
    cpV.kap1 = kapV.kap1; cpV.kap2 = kapV.kap2; % kappa_i
    cpV.kk1 = 1; cpV.kk2 = 1;    % diff coef
    cpV.lambda = model.lamV(cpV,msh.hT,CTs,pa); % penalty coef (not ghost penalty)




    %% =======================================================================
    % SOLVING U
    %=========================================================================

    %-------------------------------------------------------------------------
    % Stiffness matrix (all nodes including nodes on boundary)
    %-------------------------------------------------------------------------
    Au = getGMmPP(tris,phi,CT,msh,pa,cpU);

    %-------------------------------------------------------------------------
    % Load vector (all nodes including nodes on boundary)
    %-------------------------------------------------------------------------
    defFu = model.defFu;
    Fu = getLf(msh,pa,tris,CT,defFu);

    %-------------------------------------------------------------------------
    % BCs
    %-------------------------------------------------------------------------
    uhNX = zeros(msh.ndof,1); % column-array
    typeBC = model.bcU(); % get type of BCs
    switch typeBC
        case 1 % u=o on whole boundary
            uhNX(bNodes) = 0;
        case 2 % u=uex on whole boundary
            uhNX(bNodes) = uExNX(bNodes);
        case 3 % u=pa.bcu3 on \pt\Omg_3
            uhNX(b3Nodes) = pa.bcu3;
    end

    %----------------------------------------------------------------------
    % Solving u
    %----------------------------------------------------------------------
    Fu = Fu - Au*uhNX; % modification of F
    % LU factorization
    uhNX(iNodes) = Au(iNodes,iNodes)\Fu(iNodes); % don't care nodes on boundary
    % uhNX(iNodes) = gmres(Au(iNodes,iNodes),Fu(iNodes)); % GMRES factorization
    
    
    
    
    %% ====================================================================
    % SOLVING V
    %======================================================================
    
    %----------------------------------------------------------------------
    % Stiffness matrix (all nodes including nodes on boundary)
    %----------------------------------------------------------------------
    Av = getGMvChopp07(tris,phi,CT,msh,pa,cpV);

    %----------------------------------------------------------------------
    % Load vector (all nodes including nodes on boundary)
    %----------------------------------------------------------------------
    % in this case we use wg(u) ~ bet*u and w as bet*u, g(u) as 1
    wS = getWsep(uhNX,msh,pa.bet1,pa.bet2);
    Fv = getLvChopp07(msh,pa,tris,CT,wS);


    %----------------------------------------------------------------------
    % BCs
    %----------------------------------------------------------------------
    vhNX = zeros(msh.ndof,1); % column-array
    typeBC = model.bcV(); % get type of BCs
    switch typeBC
        case 1 % v=o on whole boundary
            vhNX(bNodes) = 0;
        case 2 % v=vex on whole boundary
            vhNX(bNodes) = vExNX(bNodes);
        case 3 % v=0 on \pt\Omg_3
            vhNX(b3Nodes) = 0;
    end

    %----------------------------------------------------------------------
    % Solving v
    %----------------------------------------------------------------------
    Fv = Fv - Av*vhNX; % modification of F
    % LU factorization
    vhNX(iNodes) = Av(iNodes,iNodes)\Fv(iNodes); % don't care nodes on boundary
    % vhNX(iNodes) = gmres(Av(iNodes,iNodes),Fv(iNodes)); % GMRES factorization
    
    
    
    
    %% ====================================================================
    % u and v in std Vh
    %======================================================================
    uhStd = interNX2STD(uhNX,msh);
    vhStd = interNX2STD(vhNX,msh);
    
    
    
    
    %% ====================================================================
    % PLOT u and v
    %======================================================================
    nf = 0; % reset every loop to be sure uh, vh plotted on the same figure
%     nf = plotNXFEM(msh,iPs,nf,'eleLabel','off','nodeLabel','off'); % only mesh
%     nf = plotNXFEM(msh,iPs,nf,uhStd,'withMesh',true,'title','uh','dim',2,'export',false); % uh
%     plotNXFEM(msh,iPs,nf,vhStd,'withMesh',true,'title','vh','dim',2,'export',false); % vh
    plotNXFEM(msh,iPs,nf,phi,'withMesh',true,'title','phi','dim',2,'export',false); % phi
    
    
%     abc = waitforbuttonpress; % wait for click

    
%     nCTs = size(iPs,3);
%     for t=1:nCTs
%         plot(iPs(1,:,t),iPs(2,:,t),'-r','LineWidth',1);
%         hold on
%     end
%     hold off
    
    pause(0.01);
    
    %% ====================================================================
    % SOLVING phi (level set function)
    % standard finite element
    %======================================================================
    
    %----------------------------------------------------------------------
    % stiffness matrix for level set
    %----------------------------------------------------------------------
%     Aphi = getGMlsChopp07(msh,pa,vhStd,dt,0.5); % use only 1 file

    % use form like in Arnold Book p.221
    del = getDells(msh,vhStd);
    Eij = getMEls(msh,pa,vhStd,del);
    Hij = getMHls(msh,pa,vhStd,del,dt,0.5);
    
    Aphi = Eij + Hij;
    
    %----------------------------------------------------------------------
    % load vector for level set
    %----------------------------------------------------------------------
    AFphi = Eij - Hij;
    phi = phi'; % row to column
    Fphi = AFphi*phi;
    
    %----------------------------------------------------------------------
    % seek phi
    %----------------------------------------------------------------------
    phi = Aphi\Fphi; % update phi
    phi = phi'; % column to row
    
%     clear msh tris CTs NCTs1 NCTs2 nodeCTs iN bN bNodes b3Nodes cpU cpV
end % for ns







