%% =======================================================================
% This file is used to find a very simple BIOFILM proposed by
% Chopp06combine
% Article: chopp06combine = 'A combined extended finite element and level set method for biofilm growth'
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
clear ; close all; % Initialization


% fixed parameters
%-------------------------------------------------------------------------
pa.degP1D = 3; % Gaussian quadrature points in 1D (polinomial functions)
pa.degP2D = 4; % Gaussian quadrature points in 2D (polinomial functions)
pa.degN = 8;    % Gaussian quadrature points in 2D (non-polynomial functions)
                % degree-#OfPoints : 1-1, 2-3, 3-4, 4-6, 5-7, 6-12,
                %                    7-13, 8-16, 9-19, 10-25, 11-27, 12-33
pa.tol = eps(1e3); % tolerance, 1e-14



JUST_TEST = 0; % only for testing on thi's machine!



%% OPTIONS
%-------------------------------------------------------------------------
model = model_chopp06combine;    % choose model. cf. file model_chopp2007.m

%% NEED TO BE CHANGED EVERY TEST CASE
savePlot = 1; % wanna save plot or not?
    testCase = '18-1'; % count the test and used to name the folder
    pathOption = 'findGood';
    moreInfo = 'TesT 18-1: find best for chopp06 ([R]restart). Slave shape test. Giong test 12 but noise<0. noise bigger, dt bigger'; % write inside file txt

%%
showPlot = 0; % wanna show plots?

% for both showPlot & savePlot
withMesh = false; 
plotGradv = 0; % plot gradient of v on cut triangles
plotContourChange = 0; % only plot the interface with time (hold on to see the track)
plotSolution = 1; % plot solution or not? (uh, vh)
    
pa.smallCut = 1;  % ignore small-support basis (1=ignore,0=no)
    pa.tH = 10; % to find the small support using (20) or (21) in arnold 2008

useFFmesh = 1; % use mesh generated by freefem++?
    reguMesh = 1; % use regular mesh or not? (only for useFFmesh=0)
    nSeg = 55;  % mesh settings (only for useFFmesh=0)
useNewton = 1; % use Newton to solve nonlinear problems?
    itol = 1e-4;
    
% ghost penalty
pa.useGP = 0; % wanna use ghost penalty term?
    pa.gam1 = 1e-6; % parameter for 1st term
    pa.gam2 = 1e-6 ; % parameter for 2nd term

% Fast marching method
useFMM = 0; % use fast marching method or not (mshdist)?
    numUseFMM = 0; % count the number of use of FMM
    alp_FMM = 0.9;
    stepUseFMM = 15; % use every 15 step (disable al_FMM method)

% SUPG
useSUPG = 1; % if 1, need to make more settings
    delEps = 1e-3;
    delSD = 0.5;

% Penalty parameters
%-------------------------------------------------------------------------
cpU.lamH = 1e8; % penalty coefficient for u (substrate)
cpV.lamH = 1e10; % penalty coefficient for v (potential)

% choose the machine to run
%-------------------------------------------------------------------------
% options: thi, gia, lehoan, blouza, gaia, google, ghost
 machine = 'google'; 
% machine = 'blouza';
% machine = 'thi';
% machine = 'ghost';
%machine = 'lehoan';


if JUST_TEST
    savePlot = 0;
    useFFmesh = 1;
    machine = 'thi';
    showPlot = 1;
    plotSolution = 1;
end



% only enable showPlot option on thi's machine
if ~strcmp(machine,'thi')
    showPlot = 0;
end



%% Model parameters
%-------------------------------------------------------------------------
pa.phiNew = 1; % using diff phi from semi circle!
if ~pa.phiNew % use semi circle
    pa.distancing = 0; % no  need to initialize
%     pa.r0 = 0.01;  % interface (like in Chopp's)
    pa.r0 = 0.05; % testing
%     pa.a = 1; % aspect ratio (p.49 Chopp 07 xfem)
else % usual in chopp06 and Carlos Conca
    pa.phiNoise = -0.03; % diff phi
    pa.phiHeight = 0.1;
    pa.distancing = 1; % make phi to be a signed distance function
end
    
pa.muS1 = 8.54932; pa.muS2 = 0;
pa.muP1 = 8.28785; pa.muP2 = 0;
%pa.bcu3 = 8.3e-6; % boundary condition for u on \pt\Omg_3
pa.bcu3 = 1e-5; % testing
cpU.kk1 = 146.88; cpU.kk2 = 183.6; % diff coef for u
cpV.kk1 = 1; cpV.kk2 = 1;    % diff coef for v
pa.f = 0.5; % volume fraction of active biomass
pa.K0 = 5e-7;

useFixedDist = 1; % use fixed distance Dirichlet condition like in Chopp?
    pa.L = 0.1; % fixed-distance of top-most Dirichlet condition
%     pa.L = 0.05; % testing

maxDay = 45; % using dt = dx/|u|
CFL = 1;



%% DOMAIN
%-------------------------------------------------------------------------
GeoDom = model.domain(); % domain


%% Mesh settings
%-------------------------------------------------------------------------
if ~useFFmesh
    disp('Mesh generated by matlab...');
    if ~reguMesh % not regular mesh?
        hEdgeMax = 2/nSeg;
        [points,edges,triangles] = initmesh(GeoDom,'hmax',hEdgeMax);    % irregular
    else
        [points,edges,triangles] = poimesh(GeoDom,nSeg,nSeg);           % regular
    end
else % using freefem nesh
    disp('Mesh generated by FreeFem++...');
   if ~pa.phiNew
        [points,edges,triangles] = ffreadmesh('./mesh/mesh_chopp06combine.msh');
    else
        [points,edges,triangles] = ffreadmesh('./mesh/mesh_chopp06_new.msh');
    end
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
phi = model.defPhi(x,y,pa); % 1 x number of points (row array)
phi(abs(phi)<pa.tol)=0; % find phi which are very small (~0) and set to 0



%% create command to run mshdist outside matlab
fprintf('Running on machine [%s]\n', machine);
switch machine
    case 'thi'
        path_nxfem = '/home/thi/Documents/nxfem/'; % thi's local machine
        path_phi = strcat(path_nxfem,'mshdist/');
        call_mshdist = strcat({'mshdist'},{' '},{path_phi},'phi'); % run in terminal
    case 'google'
        path_nxfem = '/home/thi/nxfem/';
        path_phi = strcat(path_nxfem,'mshdist/');
        call_mshdist = strcat({'mshdist'},{' '},{path_phi},'phi'); % run in terminal
    case 'ghost'
        path_nxfem = '/home/ghost/nxfem/'; 
        path_phi = strcat(path_nxfem,'mshdist/');
        call_mshdist = strcat({'mshdist'},{' '},{path_phi},'phi'); % run in terminal
    case 'gia'
        path_nxfem = '/home/gia/nxfem/'; % gia's local machine
        path_phi = strcat(path_nxfem,'mshdist/');
        call_mshdist = strcat({'mshdist'},{' '},{path_phi},'phi'); % run in terminal
    case 'lehoan'
        path_nxfem = '/home/lehoan/git/nxfem/'; % lehoan's local machine
        path_phi = strcat(path_nxfem,'mshdist/');
        call_mshdist = strcat({'mshdist'},{' '},{path_phi},'phi'); % run in terminal
    case 'blouza'
        path_nxfem = '/users/home/blouza/thi/nxfem/'; % blouza's local machine
        path_phi = strcat(path_nxfem,'mshdist/');
        call_mshdist = strcat({'/users/home/blouza/MshDist/build/mshdist'},{' '},{path_phi},'phi'); % run in terminal
    case 'gaia' % CHECK LATER!!!!
        path_nxfem = '/users/dinh/nxfem/'; % only on gaia machine
        path_phi = strcat(path_nxfem,'mshdist/');
%         call_mshdist = strcat({'mshdist'},{' '},{path_phi},'phi'); % run in terminal
end
call_mshdist = cell2mat(call_mshdist);



%% Distancing level set function (if it's not)
disp('Exporting to phi.mesh');
mshdist_w_mesh(msh,path_phi,'phi'); % export to .mesh
if pa.distancing
    tic;time=0;
    fprintf('Distancing phi... ');
    mshdist_w_sol(msh,phi,path_phi,'phi'); % export to phi.sol
    system(call_mshdist); % run 'mshdist file/to/phi' (distancing)
    phi = mshdist_r_sol(phi,path_phi,'phi'); % update phi
    fprintf('%fs\n',toc-time);
end



%% if SAVE PLOT
% Create a folder to save the plots
if savePlot
    disp('Creating folder to save plots...');
    path_machine = machine;
    if reguMesh && (~useFFmesh)
       path_regu = 'regu_';
    elseif ~useFFmesh
        path_regu = 'irregu_';
    else
        path_regu = '';
    end
    if ~useFFmesh
        path_useFF = 'matlabMesh';
    else
       path_useFF = 'FFmesh';
    end
    if useFMM 
        path_useFMM = '_wFMM'; 
    else
        path_useFMM = '_wtFMM';
    end
    if useSUPG
        path_useSUPG = '_wSUPG'; 
    else
        path_useSUPG = '_wtSUPG';
    end
    if pa.useGP
        path_useGP = '_wGP';
    else
        path_useGP = '_wtGP';
    end
    path_tris = num2str(size(triangles,2)); % number of triangles
    path_test_result = strcat(path_nxfem,'results/chopp06combine/',...
                testCase,'_',path_regu,path_useFF,path_useGP,path_useSUPG,...
                path_useFMM,'_',path_tris,'_',pathOption,'_',path_machine);
    path_test_remove = strcat({'rm -r'},{' '},{path_test_result}); % in case of duplicated folder
    path_test_remove = cell2mat(path_test_remove);
    system(path_test_remove);
    path_test_create = strcat({'mkdir'},{' '},{path_test_result}); % crfeate a new folder
    path_test_create = cell2mat(path_test_create);
    system(path_test_create);
end


if savePlot
   % Save parameters' info to file
    fileName = strcat(path_test_result,'/parameters_',num2str(nSeg),'.txt');
    fileID = fopen(fileName,'w');
        fprintf(fileID,'%s,\n',moreInfo);
        fprintf(fileID,'\n');
        fprintf(fileID,'Machine: %s,\n',machine);
        fprintf(fileID,'Model: %s,\n',model.name);
        fprintf(fileID,'hTmax: %0.10f,\n',msh.hTmax);
        fprintf(fileID,'no. triangles: %d,\n',size(triangles,2));
        fprintf(fileID,'\n');
        fprintf(fileID,'number of days: %f,\n',maxDay);
        fprintf(fileID,'__CFL: %f,\n',CFL);
        fprintf(fileID,'\n');
        fprintf(fileID,'Use FreeFem++ mesh: %d,\n',useFFmesh);
        fprintf(fileID,'Regular mesh: %d,\n',reguMesh);
        fprintf(fileID,'\n');
        fprintf(fileID,'Use small-cut: %d,\n',pa.smallCut);
        fprintf(fileID,'\n');
        fprintf(fileID,'Newton?: %d,\n',useNewton);
        fprintf(fileID,'__itol: %f,\n',itol);
        fprintf(fileID,'\n');
        fprintf(fileID,'Use FMM: %d,\n',useFMM);
        fprintf(fileID,'__al_FMM: %f,\n',alp_FMM);
        fprintf(fileID,'__numUseFMM: %d,\n',numUseFMM);
        fprintf(fileID,'__numUseFMM: %d,\n',stepUseFMM);
        fprintf(fileID,'\n');
        fprintf(fileID,'useSUPG: %d,\n',useSUPG);
        fprintf(fileID,'__delEps: %f,\n',delEps);
        fprintf(fileID,'__delSD: %f,\n',delSD);
        fprintf(fileID,'\n');
        fprintf(fileID,'use Ghost penalty: %d,\n',pa.useGP);
        fprintf(fileID,'__gam1: %f,\n',pa.gam1);
        fprintf(fileID,'__gam2: %f,\n',pa.gam2);
        fprintf(fileID,'\n');
        fprintf(fileID,'Panalty terms:\n');
        fprintf(fileID,'__lamHu: %f,\n',cpU.lamH);
        fprintf(fileID,'__lamHv: %f,\n',cpV.lamH);
        fprintf(fileID,'\n');
        fprintf(fileID,'Model parameters:\n');
        if ~pa.phiNew
            fprintf(fileID,'__r0: %f,\n',pa.r0);
        else
            fprintf(fileID,'__phiHeight: %f,\n',pa.phiHeight);
            fprintf(fileID,'__phiNoise: %f,\n',pa.phiNoise);
        end
        fprintf(fileID,'____distancing: %d,\n',pa.distancing);
        fprintf(fileID,'__muS1: %f,\n',pa.muS1);
        fprintf(fileID,'__muS2: %f,\n',pa.muS2);
        fprintf(fileID,'__muP1: %f,\n',pa.muP1);
        fprintf(fileID,'__muP2: %f,\n',pa.muP2);
        fprintf(fileID,'__bcu3: %f,\n',pa.bcu3);
        fprintf(fileID,'__cpU.k1: %f,\n',cpU.kk1);
        fprintf(fileID,'__cpU.k2: %f,\n',cpU.kk2);
        fprintf(fileID,'__cpV.k1: %f,\n',cpV.kk1);
        fprintf(fileID,'__cpV.k1: %f,\n',cpV.kk1);
        fprintf(fileID,'__f: %f,\n',pa.f);
        fprintf(fileID,'__K0: %f,\n',pa.K0);
        fprintf(fileID,'__useFixedDist: %d,\n',useFixedDist);
        fprintf(fileID,'____pa.L: %f,\n',pa.L);
%     fclose(fileID); 
end


%% =======================================================================
% EACH TIME STEP
%=========================================================================
day = 0; % count the time (days)
dt = 0; % initial, no meaning
voldSTD = zeros(msh.nStd,1); % initial vh for velocity grad v
ns=0; % number of steps


disp('Starting the loop...');
%% loop
while day < maxDay
    disp('-----------------------------');
    nf = 0; % reset every loop to be sure uh, vh plotted on the same figure
    ns = ns+1;
    fprintf('Step= %d\n',ns);
    
    %% =======================================================================
    % GET INFORMATION OF TRIANGLES
    % The same for both equations of u and v
    %=========================================================================

    % Triangles
    %-------------------------------------------------------------------------
    tris = getTriangles(phi,msh,pa); % tris has 3 factors (structure var)
    CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;

    
    % On cut triangles
    %-------------------------------------------------------------------------
    CT = getInfoCTs(CTs,phi,msh,pa); % CT has many factors (structure var)
    nodeCTs=CT.nodes; areaChildCTs=CT.areaChild; iPs=CT.iPs;

    
    % Find small-cut triangles (idx in the OLD CTs)
    %-------------------------------------------------------------------------
    if pa.smallCut
        tic;time=0;
        fprintf('Removing small cut triangles... ')
        [tris,CT] = findSmallPhi_after(msh,pa,phi,tris,CT);
%         clear CTs NCTs NCTs2 nodeCTs areaChildCTs iPs; % just in case (bad practice)
        CTs=tris.CTs;
        nodeCTs=CT.nodes; areaChildCTs=CT.areaChild; iPs=CT.iPs;
        fprintf('%fs\n',toc-time);
    end
    
    nCTs = size(iPs,3); % number of cut triangles



    %% =======================================================================
    % NODES
    %=========================================================================
    msh.nNew = nodeCTs.n; % number of new nodes (nodes around the interface)
    msh.ndof = msh.nNew + msh.nStd; % number of dofs
    msh.newNodes = getNewNodes(nodeCTs.all,msh.nStd); % vector contaning new numbering of nodes around interface, column
    msh.node = getNodes(tris,nodeCTs,msh,phi,pa); % get all nodes

    
    % boundary nodes and inner nodes
    %-------------------------------------------------------------------------
    [iN,bN] = getibNodes(msh);
    bNodes = bN.all;
    b3Nodes = bN.e3; % node on \pt\Omg_3 (top)
    
    
    
    %% top-most Dirichlet nodes
    % Only apply Dirichlet bc on the distance above top most biofilm height
    if useFixedDist
        disp('Find top-most points...');
        % nodes on interface (if exist)
        yTop = max(points(2,abs(phi)<pa.tol));
        % find the tallest
        if isempty(yTop)
            tmp(:,:) = iPs(2,:,:); % 2 x nCTs (y coor of all iPs)
            yTop = max(tmp(:));
            clear tmp;
        else
            tmp(:,:) = iPs(2,:,:); % 2 x nCTs (y coor of all iPs)
            yTop = max(yTop,max(tmp(:)));
            clear tmp;
        end
        yTop = yTop + pa.L;
        stdDNodes = find(points(2,:) >= yTop); % all nodes will be applied Dirichlet condition
        
        % include also the new nodes (msh.newNodes)
        newDNodes = msh.newNodes(intersect(nodeCTs.all,stdDNodes));
        DNodes = [stdDNodes, newDNodes']; % row array
    else
        DNodes = [];
    end
    
    
    
    %% Free Nodes!!!!!
    if ~isempty(DNodes)
        DirichletNodes = DNodes;
    else
        DirichletNodes = b3Nodes;
    end
    FreeNodes = setdiff(1:msh.ndof,DirichletNodes);
    
    
    
    %% =======================================================================
    % CONTROL PARAMETERS
    % depend on mesh and different for w and u
    % in child-functions, it's the variable 'cp'
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
    fprintf('Solving u (');
    
    % initial of iterative steps
    uold = zeros(msh.ndof,1);
    if ns ~= 1 % take the previous step as an initial (only on std nodes)
        uold(msh.node.std) = unew(msh.node.std); % initial of the iterative method
        % we only take std nodes because the size is diff at each step
    end
    % Dichlet condition
    uold(DirichletNodes,1) = pa.bcu3;
    
    difu = 100; % initial, harmless
    imax = 50; % cannot loop forever!
    step = 0;
    defFu = model.defFu; % RHS of u's equation
    defG = defGu; % nonlinear function of u (g(u))
    
    tic;time=0;
    if ~useNewton
        fprintf('Normal fixed point method)\n');
    else
        fprintf('Using Newton method)\n');
    end
    while (difu > itol) && (step<imax)
        step = step+1;
        uoldEach = getWsep(uold,msh,1,1); % analyze numSolu into each subdomain
        if ~useNewton % don't wanna use Newton method
            
            
            Au = getGMGG(tris,phi,CT,msh,pa,cpU);

            % Load vector (all nodes including nodes on boundary)
            %-------------------------------------------------------------------------
            Fu = getLfgu(msh,pa,tris,CT,uoldEach,defFu,defG.change,-pa.f*pa.muS1,-pa.f*pa.muS2);
            
            
            % BCs u
            %-------------------------------------------------------------------------
            unew = zeros(msh.ndof,1); % column-array
            unew(DirichletNodes,1) = pa.bcu3;


            % Solving u
            %----------------------------------------------------------------------
            Fu = Fu - Au*unew; % modification of F

            % LU factorization
            unew(FreeNodes) = Au(FreeNodes,FreeNodes)\Fu(FreeNodes); % don't care nodes on boundary
            % unew(FreeNodes) = gmres(Au(FreeNodes,FreeNodes),Fu(FreeNodes));   % GMRES factorization

            del = unew - uold;
        else % if use Newton
            % DF(u)del 
            coef.kk1 = pa.f*pa.muS1; coef.kk2 = pa.f*pa.muS2;
            Adel = getGMgPP(msh,pa,cpU,tris,CT,phi,uoldEach,defG.dchange,coef);
            
            % F(u)
            coef.kk1 = pa.f*pa.muS1; coef.kk2 = pa.f*pa.muS2;
            Auold = getGMgPP(msh,pa,cpU,tris,CT,phi,uoldEach,defG.changeFu,coef);
            Fdel = Auold*uold;
            
            % bc for del
            del = zeros(msh.ndof,1);
            del(DirichletNodes,1) = 0;
            
            
            % solve for del
            Fdel = Fdel - Adel*del;
            % LU
            del(FreeNodes) = Adel(FreeNodes,FreeNodes)\Fdel(FreeNodes);
%             del(FreeNodes) = gmres(Adel(FreeNodes,FreeNodes),Fdel(FreeNodes)); % GMRES
        
            % get unew
%             unew = zeros(msh.ndof,1); % reset its size every
%             unew(DirichletNodes,1) = uold(DirichletNodes,1) - del(DirichletNodes,1);
%             unew(FreeNodes,1) = uold(FreeNodes,1) - del(FreeNodes,1); % u_i+1 = u_i - del, update
            unew = uold - del;
        end
       
        delL2 = getNormL2fhNX(del,tris,CT,msh,pa);
        if ns~=1 % prevent del/0 when u0=0
            % update uold for the next step of Newton loop
            uold = unew;
            Uip1L2 = getNormL2fhNX(uold,tris,CT,msh,pa);
        else
            Uip1L2 = getNormL2fhNX(uold,tris,CT,msh,pa);
            % update uold for the next step of Newton loop
            uold = unew;
        end
        difu = delL2/Uip1L2; % |del|_L2/|u_i+1|_L2
        fprintf('___difu: %0.18f\n',difu);
    end
    fprintf('End of loop finding u...%fs\n',toc-time);
    
    
    %% ====================================================================
    % SOLVING V
    %======================================================================
    fprintf('Solving v... ');tic;time=0;
    
    % Stiffness matrix (all nodes including nodes on boundary)
    %----------------------------------------------------------------------
    Av = getGM_Chopp07v(tris,phi,CT,msh,pa,cpV);

    
    % Load vector (all nodes including nodes on boundary)
    %----------------------------------------------------------------------
    uSep = getWsep(unew,msh,1,1);
    defFv = model.defFv; % RHS of eqn of v, function handle
    Fv = getLfgu(msh,pa,tris,CT,uSep,defFv,defG.change,-pa.f*pa.muP1,-pa.f*pa.muP2);


    
    % BCs
    %----------------------------------------------------------------------
    vnew = zeros(msh.ndof,1); % column-array
    vnew(b3Nodes) = 0; % Dirichlet on \pt\Omg3
    FreeNodes = setdiff(1:msh.ndof,b3Nodes);

    
    % Solving v
    %----------------------------------------------------------------------
    Fv = Fv - Av*vnew; % modification of F
    
    % LU factorization
    vnew(FreeNodes) = Av(FreeNodes,FreeNodes)\Fv(FreeNodes); % don't care nodes on boundary
    % vnew(FreeNodes) = gmres(Av(FreeNodes,FreeNodes),Fv(FreeNodes));   % GMRES factorization
    fprintf('%fs\n',toc-time);
    
    
    
    %% ====================================================================
    % u and v in std Vh
    %======================================================================
    disp('Converting u, v to STD...');
    unewSTD = interNX2STD(unew,msh);
    vnewSTD = interNX2STD(vnew,msh);
    
    
    %% grad v
    % take at the center of triangle
    [gvnew.x,gvnew.y] = pdegrad(points,triangles,vnewSTD);
    [gvold.x,gvold.y] = pdegrad(points,triangles,voldSTD);
    
    % plot gvnew
    if plotGradv && showPlot
        disp('Show plots of grad(v)...');
    %     pdeplot(points,edges,triangles(1:3,:),'FlowData',[gvnew.x; gvnew.y]); % plot on the whole mesh

        % plot on the cut triangles only
        nf=nf+1; figure(nf);
        if withMesh
            pdemesh(points,edges,triangles); % plot mesh
            hold on;
        end
        plotInterface(msh,pa,phi,iPs); % plot the interface
        hold on
        gvx = gvnew.x; gvy = gvnew.y;
        gvx(1,NCTs1(5,:))=0; gvy(1,NCTs1(5,:))=0;
        gvx(1,NCTs2(5,:))=0; gvy(1,NCTs2(5,:))=0;
        pdeplot(points,edges,triangles(1:3,:),'FlowData',[gvx; gvy]); 
        hold off;
    end
    
    if plotGradv && savePlot
        g=figure;
        set(g, 'Visible', 'off');
        plotInterface(msh,pa,phi,iPs); % plot the interface
        hold on
        gvx = gvnew.x; gvy = gvnew.y;
        gvx(1,NCTs1(5,:))=0; gvy(1,NCTs1(5,:))=0;
        gvx(1,NCTs2(5,:))=0; gvy(1,NCTs2(5,:))=0;
        pdeplot(points,edges,triangles(1:3,:),'FlowData',[gvx; gvy]);
        titlePlot = strcat('grad vh, day = ',num2str(round(day,2)));
        title(titlePlot);
        fileName = strcat(path_test_result,'/gradv_',num2str(ns),...
                        '_','day_',num2str(round(day,2)),'.png');
        print(fileName,'-dpng','-r0');  
        hold off
        close(g); 
    end
    
    
    
     %% plot and save phi
    if showPlot
        tic;time=0;
        fprintf('Plotting phi... ');
        titlePlot = strcat('phi, day = ',num2str(round(day,2)));
        
        if plotContourChange % don't show phi's value, just show its contours with time
            if mod(ns,2)==1 % 2 step 1 time
                nf = plotNXFEM(msh,pa,phi,iPs,nf,'title',titlePlot,...
                    'withMesh',withMesh); % only mesh
            end
        else % plot phi's value also
            nf = plotNXFEM(msh,pa,phi,iPs,nf,phi,'withMesh',withMesh,...
                'title',titlePlot,'iC','b'); % phi
        end
        fprintf('%fs\n',toc-time);
    end
    
    if savePlot
        nf=nf+1; 
        f = figure(nf);
        if plotContourChange % don't show phi's value, just show its contours with time
            set(f, 'Visible', 'off');
            titlePlot = strcat('phi, day = ',num2str(round(day,2)));
            plotNXFEM(msh,pa,phi,iPs,nf,'title',titlePlot,...
                    'withMesh',withMesh,'show',false,'iC','b'); % only mesh
            hold on
        else % plot phi's value also
            set(f, 'Visible', 'off');
            titlePlot = strcat('phi, day = ',num2str(round(day,2)));
            plotNXFEM(msh,pa,phi,iPs,nf,phi,'withMesh',withMesh,...
                'title',titlePlot,'iC','b'); % phi
        end
        fileName = strcat(path_test_result,'/phi_',num2str(ns),'_',...
                        'day_',num2str(round(day,2)),'.png');

        % change size of images
        f.PaperUnits = 'inches';
        f.PaperPosition = [0 0 8 6];
        print(fileName,'-dpng','-r0');
        if ~plotContourChange
            close(f);
        end 
    end
    
    
    
    if dt > 15
        break;
    end
    
    
    
    %% ====================================================================
    % PLOT u and v
    %======================================================================
    if showPlot && plotSolution
%         nf = plotNXFEM(msh,pa,phi,iPs,nf,'eleLabel','off','nodeLabel','off'); % only mesh

        disp('Plotting u...');
        titlePlot = strcat('uh, day = ',num2str(round(day,2)));
        nf = plotNXFEM(msh,pa,phi,iPs,nf,unewSTD,'withMesh',withMesh,...
                    'title',titlePlot,'iC','b'); % uh

        disp('Plotting v...');
        titlePlot = strcat('vh, day = ',num2str(round(day,2)));
        nf = plotNXFEM(msh,pa,phi,iPs,nf,vnewSTD,'withMesh',withMesh,...
                    'title',titlePlot,'iC','b'); % vh
    end % end if showPlot
    
    
    
    %% SAVE PLOT
    if savePlot
        tic;time=0;
        fprintf('Saving plots... ');
        
        if plotSolution
            % uh
            nf=nf+1;
            g=figure(nf);
            set(g, 'Visible', 'off');
            titlePlot = strcat('uh, day = ',num2str(round(day,2)));
            plotNXFEM(msh,pa,phi,iPs,nf,unewSTD,'withMesh',withMesh,...
                            'title',titlePlot,'show',false,'iC','b'); % uh
            fileName = strcat(path_test_result,'/uh_',num2str(ns),...
                            '_','day_',num2str(round(day,2)),'.png');
            % change size of images
            g.PaperUnits = 'inches';
            g.PaperPosition = [0 0 8 6];
            print(fileName,'-dpng','-r0');
            close(g);

            % vh
            nf=nf+1;
            g=figure(nf);
            set(g, 'Visible', 'off');
            titlePlot = strcat('vh, day = ',num2str(round(day,2)));
            plotNXFEM(msh,pa,phi,iPs,nf,vnewSTD,'withMesh',withMesh,...
                            'title',titlePlot,'show',false,'iC','b'); % vh
            fileName = strcat(path_test_result,'/vh_',num2str(ns),...
                            '_','day_',num2str(round(day,2)),'.png');
            % change size of images
            g.PaperUnits = 'inches';
            g.PaperPosition = [0 0 8 6];
            print(fileName,'-dpng','-r0');
            close(g);
        end
        
        fprintf('%fs\n',toc-time);
    end
    
    
    
    %% ====================================================================
    % SOLVING phi (level set function)
    % standard finite element
    %======================================================================
    disp('Solving level set phi...');
    
    % get del_T
    if useSUPG
        delOld = getDellsT(msh,gvold,delEps,delSD); % Arnold's book p.223
        delNew = getDellsT(msh,gvnew,delEps,delSD);
    else
        delOld = zeros(1,size(msh.t,2)); % without SUPG
        delNew = delOld;
    end
    
    
    % dt
    maxGradV = max(abs(gvnew.x) + abs(gvnew.y));
    dt = CFL*msh.hTmax/maxGradV;
    day = day+dt;
    fprintf('day: %f\n',day);
    fprintf('dt: %f\n',dt);
    
    
    
    if savePlot && (ns==1)
        fprintf(fileID,'The first dt: %f,\n',dt);
        fprintf(fileID,'\n');
    end
    if ns==1
        fprintf('The first dt: %f\n',dt);
    end
    
    
    % stiffness matrix for level set
    %----------------------------------------------------------------------
    tic;time=0;
    fprintf('Get stiffness matrix Enew, Hnew... ');
    Enew = getMElsGP(msh,pa,gvnew,delNew,1);
    Hnew = getMHlsGP(msh,pa,gvnew,delNew,dt*0.5);
    mI = speye(msh.nStd); % identity matrix
    Aphi = mI + Enew^(-1)*Hnew;
    fprintf('%fs\n',toc-time);
    
    
    % load vector for level set
    %----------------------------------------------------------------------
    tic;time=0;
    fprintf('Get load vector (Eold, Hold, Afphi)... ');
    Eold = getMElsGP(msh,pa,gvold,delOld,1);
    Hold = getMHlsGP(msh,pa,gvold,delOld,dt*0.5);
    AFphi = mI - Eold^(-1)*Hold;
    phi = phi'; % row to column
    Fphi = AFphi*phi;
    fprintf('%fs\n',toc-time);
    
    
    % seek phi
    %----------------------------------------------------------------------
    disp('Updating phi...');
    phi = Aphi\Fphi; % update phi
    phi = phi'; % column to row
    
    
    
    
    %% update v
    disp('Updating v...');
    voldSTD = vnewSTD;
    
    
    
    
    %% Reinitialization
    %----------------------------------------------------------------------
    norm_gradphi = getNormL2GfhSTD(msh,phi); % ||gradPhi||_L2
    fprintf('|1-norm_gradphi| = %f\n',abs(1-norm_gradphi));
    
%     if useFMM && abs(1-norm_gradphi) > alp_FMM && numUse <=1
     if useFMM && abs(1-norm_gradphi) > alp_FMM
%    if useFMM && (mod(ns,stepUseFMM)==0) % every stepUseFMM step
        disp('Starting to use FMM...');
        mshdist_w_sol(msh,phi,path_phi,'phi'); % export to phi.sol
        system(call_mshdist); % run 'mshdist file/to/phi' (redistancing)
        phi = mshdist_r_sol(phi,path_phi,'phi'); % update phi
        numUseFMM = numUseFMM + 1;
    end
    
end % for ns


%% save info file
if savePlot
   % Save parameters' info to file
%     fileName = strcat(path_test_result,...
%         '/parameters_',num2str(size(triangles,2)),'.txt');
%     fileID = fopen(fileName,'w');
        fprintf(fileID,'\n');
        fprintf(fileID,'the last dt: %f,\n',dt);
        fprintf(fileID,'\n');
    fclose(fileID); 
end


%% CLEAN UP
% Fix conflict in git
system('rm -r mshdist/phi.mesh');
close all; % close all figures in this test



function delT = getDellsT(msh,gP,eps,SD)
    % This function 'guest' formula presented in Arnold's book, page 223 (and also cf. (7.14))
    % If in future, we need it more than once, I will put it in a separated file
    % Input: - msh.hT : h on each triangle : 1 x nTs
    %        - gP: gP.x (1 x nTs), gP.y (1 x nTs)
    %        - a control number SD \in O(1) (cf. (7.14))
    %        - given small eps > 0 (cf. 7.14), at page 223, he took eps=1e-3
    % Output: a vector 1 x nTs
    
    maxGradV = max(abs(gP.x) + abs(gP.y));
    delT = SD*msh.hT/max(eps,maxGradV);
end