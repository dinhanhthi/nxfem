%% =======================================================================
% This file is used to find a very simple BIOFILM proposed by
% Chopp06combine
% Article: chopp06combine = "A combined extended finite element and level set method for biofilm growth"
% This is a SYSTEM of LINEAR equations u,v (u=substrate=S, v=potential=Phi)
% ------------------------------------------------------------------------
% PURPOSE: verifying the results given in chopp06combine
% Related files: 
%   - model_chopp06combine.m: contains all info for the model
% ------------------------------------------------------------------------
% RESULT: 
% ------------------------------------------------------------------------
% Last modified:
% ------------------------------------------------------------------------

%% =======================================================================
% DOMAIN: [0,0.5]X[0,0.1], inital interface: half of circle ((0.25,0);r0)
% all meters are in mm
%-------------------------------------------------------------------------
% MODELS:
% r=sqrt((x-0.25)^2+y^2)-r0 the interface
% Omg1 = Omg_b, Omg2 = Omg_f
%-------------------------------------------------------------------------
% Equation of u (substrate S):
% -nabla(Ds*nabla(u)) + f*muS*u/(K0+u) = 0  in Omega
% [u]=[Ds*nabla_n(u)]=0     on Gamma
% u = 8.3e-6   on \partial\Omega_3
% nabla_n u = 0 elsewhere
%-------------------------------------------------------------------------
% Equation of v (potential Phi):
% -nabla(nabla(v)) = -f*muP*u/(K0+u) in Omega
% v=nabla_n(v)=0 on Gamma 
% v=0 on \partial\Omega_3 (checked again!)
% nabla_n v = 0 elsewhere
%-------------------------------------------------------------------------
% PARAMETERS:
% Ds = 146.88 (Omg1), 183.6 (Omg2) (DSb=0.8DSf)
% muS = 8.54932 (Omg1), 0 (Omg2)
% muP = 8.28785 (Omg1), 0 (Omg2)
%=========================================================================



%% =======================================================================
% PARAMETERS
% Note that parameters are different for equations of u and v
%=========================================================================

addpath(genpath('func'));   % add all necessary functions

% fixed parameters
%-------------------------------------------------------------------------
pa.degP1D = 3;              % Gaussian quadrature points in 1D (polinomial functions)
pa.degP2D = 4;              % Gaussian quadrature points in 2D (polinomial functions)
pa.degN = 8;    % Gaussian quadrature points in 2D (non-polynomial functions)
                % degree-#OfPoints : 1-1, 2-3, 3-4, 4-6, 5-7, 6-12,
                %                    7-13, 8-16, 9-19, 10-25, 11-27, 12-33
pa.tol = eps(1e3);          % tolerance, 1e-14


% MODEL
%-------------------------------------------------------------------------
model = model_chopp06combine;    % choose model. cf. file model_chopp2007.m
showPlot = 1;
    withMesh = true;

pa.smallCut = 0;  % ignore small-support basis (1=ignore,0=no)
pa.tH = 10; % to find the small support using (20) or (21) in arnold 2008

pa.reguMesh = 0;            % use regular mesh or not?
nSeg = 127;  % mesh settings
useNewton = 0; % use Newton to solve nonlinear problems?

% Ghost penalty
%-------------------------------------------------------------------------
pa.useGP = 1;   % wanna use ghost penalty term?
    pa.gam1 = 1e-6;             % parameter for 1st term
    pa.gam2 = 1e-6 ;            % parameter for 2nd term


% Model parameters
%-------------------------------------------------------------------------
% pa.r0 = 0.01;  % interface
pa.r0 = 0.05; % testing
    pa.distancing = 0; % make phi to be a signed distance function
pa.muS1 = 8.54932; pa.muS2 = 0;
pa.muP1 = 8.28785; pa.muP2 = 0;
pa.bcu3 = 8.3e-6; % boundary condition for u on \pt\Omg_3
cpU.kk1 = 146.88; cpU.kk2 = 183.6; % diff coef for u
cpV.kk1 = 1; cpV.kk2 = 1;    % diff coef for v
pa.f = 0.5; % volume fraction of active biomass
pa.K0 = 5e-7;


% Penalty parameters
%-------------------------------------------------------------------------
cpU.lamH = 1e5; % penalty coefficient for u (substrate)
cpV.lamH = 1e5; % penalty coefficient for v (potential)



% DOMAIN
%-------------------------------------------------------------------------
GeoDom = model.domain(); % domain


% Mesh settings
%-------------------------------------------------------------------------
if ~pa.reguMesh % not regular mesh?
    hEdgeMax = 2/nSeg;
    [points,edges,triangles] = initmesh(GeoDom,'hmax',hEdgeMax);    % irregular
else
    [points,edges,triangles] = poimesh(GeoDom,nSeg,nSeg);           % regular
end

msh.p = points; msh.t = triangles; msh.e = edges;   % save to msh
x = points(1,:);    % x-coordinate of points
y = points(2,:);    % y-coordinate of points

% diameter (longest side) of each triangle: 1 x nTs
msh.hT = getDiam(msh);              % 1 x number of triangles
msh.hTmax = max(msh.hT);            % maximum of all diameters
msh.nStd = size(points,2);          % number of standard nodes


% Level set function (INITIAL)
%-------------------------------------------------------------------------
phi = model.defPhi(x,y,pa);         % 1 x number of points (row array)
% phi(abs(phi)<pa.tol)=0;             % find phi which are very small (~0) and set to 0



% level set settings
%-------------------------------------------------------------------------
useFMM = 1; % use fast marching method or not (mshdist)?
    numUse = 0; % count the number of use of FMM
    alp_FMM = 0.1;
useSUPG = 1; % if 1, need to make more settings
    delEps = 1e-3;
    delSD = 0.5;
    
path_nxfem = '/home/thi/Dropbox/git/nxfem/'; % thi's local machine
% path_nxfem = '/users/dinh/nxfem/'; % only on gaia machine
path_phi = strcat(path_nxfem,'mshdist/');
call_mshdist = strcat({'mshdist'},{' '},{path_phi},'phi'); % run in terminal
call_mshdist = cell2mat(call_mshdist);


%% Distancing level set function (if it's not)
disp("Exporting to phi.mesh");
mshdist_w_mesh(msh,path_phi,'phi'); % export to .mesh
if pa.distancing
    disp("Distancing phi...");
    mshdist_w_sol(msh,phi,path_phi,'phi'); % export to phi.sol
    system(call_mshdist); % run 'mshdist file/to/phi' (distancing)
    phi = mshdist_r_sol(phi,path_phi,'phi'); % update phi
end



%% =======================================================================
% EACH TIME STEP
%=========================================================================
% pA = model.pa;
% pA = pA();
maxStep = 50; % using dt = dx/|u|
vhOld = zeros(msh.nStd,1); % initial vh for velocity grad v

disp("Starting the loop...");
%% loop
for ns = 1:maxStep
    disp("-----------------------------");
    Xdisp = ['step = ', num2str(ns)];
    disp(Xdisp);
    
    %% =======================================================================
    % GET INFORMATION OF TRIANGLES
    % The same for both equations of u and v
    %=========================================================================

    % Triangles
    %-------------------------------------------------------------------------
    tris = getTriangles(phi,msh,pa);    % tris has 3 factors (structure var)
    CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;

    
    % On cut triangles
    %-------------------------------------------------------------------------
    CT = getInfoCTs(CTs,phi,msh,pa);    % CT has many factors (structure var)
    nodeCTs=CT.nodes; areaChildCTs=CT.areaChild; iPs=CT.iPs;

    
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
    msh.nNew = nodeCTs.n;   % number of new nodes (nodes around the interface)
    msh.ndof = msh.nNew + msh.nStd;         % number of dofs
    msh.newNodes = getNewNodes(nodeCTs.all,msh.nStd);       % vector contaning new numbering of nodes around interface, column
    msh.node = getNodes(tris,nodeCTs,msh,phi,pa);           % get all nodes

    
    % boundary nodes and inner nodes
    %-------------------------------------------------------------------------
    [iN,bN] = getibNodes(msh);
    bNodes = bN.all; iNodes=iN.all;
    b3Nodes = bN.e3;        % node on \pt\Omg_3 (top)
    
    
    disp("Plotting phi...");
    %% plot phi
    nf = 0; % reset every loop to be sure uh, vh plotted on the same figure
%     abc = waitforbuttonpress; % wait for click
    titlePlot = strcat('phi, step = ',num2str(ns));
    if showPlot
        nf = plotNXFEM(msh,iPs,nf,phi,'withMesh',withMesh,'title',titlePlot,'dim',2,'export',false); % phi
    end
    nCTs = size(iPs,3);
    hold on
    for t=1:nCTs
        plot(iPs(1,:,t),iPs(2,:,t),'-r','LineWidth',1);
        hold on
    end
    hold off
%     pause(0.01);
    
    
    %% =======================================================================
    % CONTROL PARAMETERS
    % depend on mesh and different for w and u
    % in child-functions, it's the variable "cp"
    %=========================================================================
    hTCTs = msh.hT(CTs(5,:));
    
    kapU = model.kapU(cpU,CT,pa);
    cpU.kap1 = kapU.kap1; cpU.kap2 = kapU.kap2; % kappa_i
    cpU.lambda = model.lamU(cpU,hTCTs,CT,pa); % penalty coef (not ghost penalty)

    kapV = model.kapV(cpV,CT,pa);
    cpV.kap1 = kapV.kap1; cpV.kap2 = kapV.kap2; % kappa_i
    cpV.lambda = model.lamV(cpV,hTCTs,CT,pa); % penalty coef (not ghost penalty)




    %% =======================================================================
    % SOLVING U
    %=========================================================================

    
    % Stiffness matrix (all nodes including nodes on boundary)
    %-------------------------------------------------------------------------
    disp("Solving u...");
    
    % initial of iterative steps
    if ns==1
        uold = zeros(msh.ndof,1);
    end
    itol = 1e-2;
    difu = 100; % initial, useless
    imax = 50;
    step = 0;
    typeBC = model.bcU();   % get type of BCs
    defFu = model.defFu;
    defG = defGu;
    
    if ~useNewton
        disp("Normal fixed point method.");
    else
        disp("Using Newton method");
    end
    while (difu > itol) && (step<imax)
        step = step+1;
        uoldEach = getWsep(uold,msh,1,1); % analyze numSolu into each subdomain
        if ~useNewton % don't wanna use Newton method
            Au = getGMGG(tris,phi,CT,msh,pa,cpU);

            % Load vector (all nodes including nodes on boundary)
            %-------------------------------------------------------------------------
            Fu = getLfgu(msh,pa,tris,CT,uoldEach,defFu,defG.change,-pa.f*pa.muS1,-pa.f*pa.muS2);



            % BCs
            %-------------------------------------------------------------------------
            unew = zeros(msh.ndof,1); % column-array
            switch typeBC
                case 1 % u=o on whole boundary
                    unew(bNodes) = 0;
                case 2  % u=uex on whole boundary
                    unew(bNodes) = uExNX(bNodes);
                case 3 % u=pa.bcu3 on \pt\Omg_3
                    unew(b3Nodes) = pa.bcu3;
            end


            % Solving u
            %----------------------------------------------------------------------
            Fu = Fu - Au*unew; % modification of F

            % LU factorization
            unew(iNodes) = Au(iNodes,iNodes)\Fu(iNodes); % don't care nodes on boundary
            % unew(iNodes) = gmres(Au(iNodes,iNodes),Fu(iNodes));   % GMRES factorization

            del = unew - uold;
            uold = unew; % update for the next step
        end
        delL2 = getNormL2nxfem(del,tris,CT,msh,pa);
        Uip1L2 = getNormL2nxfem(uold,tris,CT,msh,pa);
        difu = delL2/Uip1L2; % |del|_L2/|u_i+1|_L2
%         difu = delL2;
        fprintf('difu: %0.18f\n',difu);
    end
    disp("End of loop finding u");
    
    
    %% ====================================================================
    % SOLVING V
    %======================================================================
    disp("Solving v...");
    
    % Stiffness matrix (all nodes including nodes on boundary)
    %----------------------------------------------------------------------
    Av = getGM_Chopp07v(tris,phi,CT,msh,pa,cpV);

    
    % Load vector (all nodes including nodes on boundary)
    %----------------------------------------------------------------------
    uSep = getWsep(unew,msh,1,1);
    Fv = getLfgu(msh,pa,tris,CT,uSep,defFu,defG.change,-pa.f*pa.muP1,-pa.f*pa.muP2);


    
    % BCs
    %----------------------------------------------------------------------
    vnew = zeros(msh.ndof,1);           % column-array
    typeBC = model.bcV();               % get type of BCs
    switch typeBC
        case 1          % v=o on whole boundary
            vnew(bNodes) = 0;
        case 2          % v=vex on whole boundary
            vnew(bNodes) = vExNX(bNodes);
        case 3          % v=0 on \pt\Omg_3
            vnew(b3Nodes) = 0;
    end

    
    % Solving v
    %----------------------------------------------------------------------
    Fv = Fv - Av*vnew;      % modification of F
    
    % LU factorization
    vnew(iNodes) = Av(iNodes,iNodes)\Fv(iNodes);            % don't care nodes on boundary
    
    % vnew(iNodes) = gmres(Av(iNodes,iNodes),Fv(iNodes));   % GMRES factorization
    
    
    
    
    %% ====================================================================
    % u and v in std Vh
    %======================================================================
    disp("Converting u, v to STD...");
    uhSTD = interNX2STD(unew,msh);
    vhSTD = interNX2STD(vnew,msh);
    
    
    
    
    %% ====================================================================
    % PLOT u and v
    %======================================================================
%     nf = 0; % reset every loop to be sure uh, vh plotted on the same figure
%     nf = plotNXFEM(msh,iPs,nf,'eleLabel','off','nodeLabel','off'); % only mesh
    
    titlePlot = strcat('uh, step = ',num2str(ns));
    nf = plotNXFEM(msh,iPs,nf,uhSTD,'withMesh',true,'title',titlePlot,'dim',2,'export',false); % uh
    
    titlePlot = strcat('vh, step = ',num2str(ns));
    nf = plotNXFEM(msh,iPs,nf,vhSTD,'withMesh',true,'title',titlePlot,'dim',2,'export',false); % vh

%     plotNXFEM(msh,iPs,nf,phi,'withMesh',true,'title','phi','dim',2,'export',false); % phi
  
    hold on
    for t=1:nCTs
        plot(iPs(1,:,t),iPs(2,:,t),'-r','LineWidth',1);
        hold on
    end
    hold off

%     abc = waitforbuttonpress; % wait for click
%     pause(0.01);
    
    %% ====================================================================
    % SOLVING phi (level set function)
    % standard finite element
    %======================================================================
    disp("Solving level set phi...");
    
    % get del_T
    if useSUPG
        delOld = getDellsT(msh,vhOld,delEps,delSD); % Arnold's book p.223
        delNew = getDellsT(msh,vhSTD,delEps,delSD);
    else
        delOld = zeros(1,size(msh.t,2)); % without SUPG
        delNew = delOld;
    end
    
    
    % dt
    dt = msh.hTmax/max(delNew);
    
    
    % stiffness matrix for level set
    %----------------------------------------------------------------------
    Enew = getMEls_gP(msh,pa,vhSTD,delNew,1);
    Hnew = getMHls_gP(msh,pa,vhSTD,delNew,dt*0.5);
    mI = speye(msh.nStd); % identity matrix
    Aphi = mI + Enew^(-1)*Hnew;
    
    
    % load vector for level set
    %----------------------------------------------------------------------
    Eold = getMEls_gP(msh,pa,vhOld,delOld,1);
    Hold = getMHls_gP(msh,pa,vhOld,delOld,dt*0.5);
    AFphi = mI - Eold^(-1)*Hold;
    phi = phi';             % row to column
    Fphi = AFphi*phi;
    
    
    % seek phi
    %----------------------------------------------------------------------
    disp("Updating phi...");
    phi = Aphi\Fphi;        % update phi
    phi = phi';             % column to row
    
    
    %% update v
    disp("Updating v...");
    vhOld = vhSTD;
    
    
    %% Reinitialization
    %----------------------------------------------------------------------
    norm_gradphi = getNormL2Gstd(msh,phi); % ||gradPhi||_L2
    todisplayed = ['|1-norm_gradphi| = ',num2str(abs(1-norm_gradphi))];
    disp(todisplayed);
    
%     if useFMM && abs(1-norm_gradphi) > alp_FMM && numUse <=1
    if useFMM && abs(1-norm_gradphi) > alp_FMM
        disp("Starting to use FMM...");
        mshdist_w_sol(msh,phi,path_phi,'phi'); % export to phi.sol
        system(call_mshdist); % run 'mshdist file/to/phi' (redistancing)
        phi = mshdist_r_sol(phi,path_phi,'phi'); % update phi
        numUse = numUse + 1;
    end
    
    
    
end % for ns


function delT = getDellsT(msh,vold,eps,SD)
    % find ||grad v||_inf on each triangle
    % This function "guest" formula presented in Arnold's book, page 223 (and also cf. (7.14))
    % If in future, we need it more than once, I will put it in a separated file
    % old file: getDells.m (output delT is a scalar)
    % Input: - msh.hT : h on each triangle : 1 x nTs
    %        - vold to find grad v
    %        - a control number SD \in O(1) (cf. (7.14))
    %        - given small eps > 0 (cf. 7.14), at page 223, he took eps=1e-3
    % Output: a vector 1 x nTs
    
    tris = msh.t;
    nTs = size(tris,2);
    delT = zeros(1,nTs);
    gP = getGradPhi(tris,msh); % 2 coordinates x 3 vertices x nTris
    for t=1:nTs
        v1t = abs(vold(tris(1,t)))*( abs(gP(1,1,t)) + abs(gP(2,1,t)) ); % vertex 1
        v2t = abs(vold(tris(2,t)))*( abs(gP(1,2,t)) + abs(gP(2,2,t)) ); % vertex 2
        v3t = abs(vold(tris(3,t)))*( abs(gP(1,3,t)) + abs(gP(2,3,t)) ); % vertex 3
        delT(1,t) = max([eps, v1t, v2t, v3t]);
    end
    
    delT = SD*msh.hT ./ delT;
end