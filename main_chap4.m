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
% RESULT DIFFERENT FROM THE PAST, NEED TO CHECK IT AGAIN LATER!!!!!
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
clear ; close all; % Initialization



%% Fixed parameters
pa.degP1D = 3; % Gaussian quadrature points in 1D (for polinomial functions)
pa.degP2D = 4; % Gaussian quadrature points in 2D (for polinomial functions)
pa.degN = 8; % Gaussian quadrature points in 2D (for non-polynomial functions)
% degree-#OfPoints : 1-1, 2-3, 3-4, 4-6, 5-7, 6-12, 7-13, 8-16, 9-19, 10-25, 11-27, 12-33
pa.tol = eps(1e3); % tolerance, 1e-14



%% Model
model = model_article1;
GeoDom = model.domain(); % domain
defG = defGu; % def of gu, cf. defGu.m


%% SETTINGS
findCR = 0; % wanna find the CR?
    %     numSegCR = [16, 32, 64, 128]; % only works with findCR=1
    numSegCR = [36, 56, 86, 126];
    numSegPlot = 37; % only for plotting, findCR=0
savePlot = 0; % 1 = export figures to files (and without plotting them)
    testCase = '1';
    pathOption = '';
    moreInfo = 'Test 1: First test'; % write inside file txt
showPlot = 1; % wanna plot or not the solution?
    nf = 0; % counter of figures (plot each plot in a separated figure)
    withMesh = false; % plot with mesh?
    
pa.smallCut = 0; % ignore small-support basis (1=ignore,0=no)
    pa.tH = 1e2; % to find the small support using (20) or (21) in arnold 2008
    
reguMesh = 1; % regular or irregular mesh?

useNewton = 0; % use Newton method for solving v?
    itol = 1e-3;
    imax = 50; % number of iterative

% Penalty (goes with \int [][])
cpW.lamH = 1e8; % penalty coefficient for w
cpV.lamH = 1e8; % penalty coefficient for v

% ghost penalty
pa.useGP = 1; % wanna use ghost penalty term?
    pa.gam1 = 1e-6; % parameter for 1st term
    pa.gam2 = 1e-6 ; % parameter for 2nd term



%% Model's parameters
cpW.kk1 = 1; cpW.kk2 = 100;
cpV.kk1 = 0.5; cpV.kk2 = 100;
% just for findDef?ex
pa.alp1 = cpW.kk1; pa.alp2 = cpW.kk2;
pa.bet1 = cpV.kk2; pa.bet2 = cpV.kk2;

% (cpV.kk2 doesn't take affect because v=0 in Omg2)
pa.r0 = 0.6; % interface
pa.lamSys = 1; % coef lam in system settings



%% choose the machine to run
% options: thi, gia, lehoan, blouza, gaia, google, ghost
% machine = 'google'; 
% machine = 'blouza';
machine = 'thi';
% machine = 'ghost';
% machine = 'lehoan';

% only enable showPlot option on thi's machine
if ~strcmp(machine,'thi')
    showPlot = 0;
end



%% create command to run mshdist outside matlab
fprintf('Running on machine [%s]\n', machine);
switch machine
    case 'thi'
        path_nxfem = '/home/thi/Documents/nxfem/'; % thi's local machine
    case 'google'
        path_nxfem = '/home/thi/nxfem/';
    case 'ghost'
        path_nxfem = '/home/ghost/nxfem/'; 
    case 'gia'
        path_nxfem = '/home/gia/nxfem/'; % gia's local machine
    case 'lehoan'
        path_nxfem = '/home/lehoan/git/nxfem/'; % lehoan's local machine
    case 'blouza'
        path_nxfem = '/users/home/blouza/thi/nxfem/'; % blouza's local machine
    case 'gaia' % CHECK LATER!!!!
        path_nxfem = '/users/dinh/nxfem/'; % only on gaia machine
end



%% Dependent parameters
if findCR == 1
    numSeg = numSegCR;
    disp('Find the convergence rate...');
else
    numSeg = numSegPlot;
    disp('Only ONE step...');
end
nStep = size(numSeg,2);



%% if SAVE PLOT
% Create a folder to save the plots
if savePlot
    disp("Creating folder to save plots...");
    path_machine = machine;
    if reguMesh && (~useFFmesh)
       path_regu = 'regu_';
    elseif ~useFFmesh
        path_regu = 'irregu_';
    else
        path_regu = '';
    end
    path_test_result = strcat(path_nxfem,'results/chap4/',...
                testCase,'_',path_regu,'_',pathOption,'_',path_machine);
    path_test_remove = strcat({'rm -r'},{' '},{path_test_result}); % in case of duplicated folder
    path_test_remove = cell2mat(path_test_remove);
    system(path_test_remove);
    path_test_create = strcat({'mkdir'},{' '},{path_test_result}); % crfeate a new folder
    path_test_create = cell2mat(path_test_create);
    system(path_test_create);
end



%% For finding CR
hTarray = zeros(1,nStep);
nDOFsArray = zeros(1,nStep);
wL2array = zeros(1,nStep);
vL2array = zeros(1,nStep);
uL2array = zeros(1,nStep);



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
    
    
    %% Exact solution in stdFEM
    disp('Exact solutions...');
    % w
    defExSolw = model.defWex; % function handle
    wExSTD = exInStd(defExSolw,msh,pa);
    % u
    defExSolu = model.defUex; % function handle
    uExSTD = exInStd(defExSolu,msh,pa);
    % v
    defExSolv = model.defVex; % function handle
    vExSTD = exInStd(defExSolv,msh,pa);
    
    
    %% Exact solution in NXFEM
    wExNX = interSTD2NX(wExSTD,msh); % column array
    % u
    uExNX = interSTD2NX(uExSTD,msh); % column array
    % v
    vExNX = interSTD2NX(vExSTD,msh); % column array
    
    
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
    Fw = getLf(msh,pa,tris,CT,defFw);
    
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
    vold = zeros(msh.ndof,1); % initial solution of v
    difu = 100; % initial, harmless
    step = 0;
    defFvGu = model.defFv; % (x,y,pa,sub,defG)
    defFv = @(xx,yy,pa,sub) defFvGu(xx,yy,pa,sub,defG.change);
    typeBC = model.bcV(); % get type of BCs
    wEach = getWsep(whNX,msh,1,1); % w in each subdomain
    
    if ~useNewton
        fprintf('Normal fixed point method)\n');
    else
        fprintf('Using Newton method)\n');
    end
    
    while (difu > itol) && (step<imax)
        step = step+1;
        
        % bet/(alp*lam)v
        voldEach = getWsep(vold,msh,pa.bet1/(pa.alp1*pa.lamSys),...
                            pa.bet2/(pa.alp2*pa.lamSys));
                        
        % w-bet/(alp*lam)v
        wvEach.omg1 = wEach.omg1 - voldEach.omg1;
        wvEach.omg2 = wEach.omg2 - voldEach.omg2;
        wvEach.ct1 = wEach.ct1 - voldEach.ct1;
        wvEach.ct2 = wEach.ct2 - voldEach.ct2;
                        
        if ~useNewton % don't use Newton
            
            % -grad*grad - lam*g(wold-bet/(alp*lam)vold)*v*phi
            coef.kk1 = pa.lamSys; coef.kk2 = pa.lamSys;
            Av = getGMvAA(msh,pa,tris,CT,phi,wvEach,cpV,defG.change);
                       
            Fv = getLf(msh,pa,tris,CT,defFv);
            
            vnew = zeros(msh.ndof,1); % zero initial uh for each step
            switch typeBC
                case 1 % u=o on whole boundary
                    vnew(bNodes) = 0;
                case 2 % u=uex on whole boundary
                    vnew(bNodes) = vExNX(bNodes);
            end
            
            Fv = Fv - Av*vnew; % modification of F
            vnew(iNodes) = Av(iNodes,iNodes)\Fv(iNodes);
            % vnew(iNodes) = gmres(Av(iNodes,iNodes),Fv(iNodes)); 
            
            del = vnew - vold;
            vold = vnew;
        else
            % DF(u)
            Adel = getGMvAANewton(msh,pa,tris,CT,phi,wvEach,voldEach,cpV,defG);
           
            
            % F(u)
            Av = getGMvAA(msh,pa,tris,CT,phi,wvEach,cpV,defG.change); % like Av in normal iterative method
            getGMvAA(tris,phi,voldEach,wEach,CT,msh,pa,cpV,cpW); 
            Fv = getLf(msh,pa,tris,CT,defFv); % like Fu in normal iterative method
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
                
    
    
    %% Solution NX in STD (for plotting)
    whSTD = interNX2STD(whNX,msh);
    uhSTD = interNX2STD(uhNX,msh);
    vhSTD = interNX2STD(vhNX,msh);
                
                
                
    %% ERRORS
    fprintf('Get errors...');
    wL2 = getNormL2uhuexNX(msh,pa,tris,CT,uhNX,defExSolw,defPhi);
    vL2 = getNormL2uhuexNX(msh,pa,tris,CT,vhNX,defExSolv,defPhi);
    uL2 = getNormL2uhuexNX(msh,pa,tris,CT,uhNX,defExSolu,defPhi);
    
    
    
    %% For finding CR
    hTarray(z) = msh.hTmax;
    nDOFsArray(z) = msh.ndof; 
    wL2array(z) = wL2;
    vL2array(z) = vL2;
    uL2array(z) = uL2;
end



%% Display CR or plotting
if nStep>1
    %% find CR
    disp('Find the CR...');
    tmp = polyfit(log(hTarray),log(wL2array),1);
    orderL2w = tmp(1);
    
    tmp = polyfit(log(hTarray),log(vL2array),1);
    orderL2v = tmp(1);
    
    tmp = polyfit(log(hTarray),log(uL2array),1);
    orderL2u = tmp(1);
    
    
    %% Save figure of CR
    if savePlot
        disp('Saving plot of CR...');
        
        %% Save figure of CR
        nf = nf+1;
        f = figure(nf);
        set(f, 'Visible', 'off');
        plot(log(hTarray),log(wL2array),'-r.',...
            log(hTarray),log(uL2array),'-b.',...
            log(hTarray),log(vL2array),'-m.');
        legend('w','u','v');
        xlabel('log(h)'); 
        ylabel('log(L2 error)');
        fileName = strcat(path_test_result,'/CR_',testCase,'.png');
        % change size of images
        f.PaperUnits = 'inches';
        f.PaperPosition = [0 0 8 6];
        print(fileName,'-dpng','-r0');
        if ~plotContourChange
            close(f);
        end
        
        
        %% save errors as file
        % see file results\main_nxfem\errors.txt
        cr = zeros(3,nStep);
        for i=2:nStep
           cr(1,i) = log(wL2array(i)/wL2array(i-1))/log(hTarray(i)/hTarray(i-1)); % w
           cr(2,i) = log(uL2array(i)/uL2array(i-1))/log(hTarray(i)/hTarray(i-1)); % u
           cr(3,i) = log(vL2array(i)/vL2array(i-1))/log(hTarray(i)/hTarray(i-1)); % v
        end
        errFileMat = [hTarray;wL2array;cr(1,:);uL2array;cr(2,:);vL2array;cr(3,:)];
        fileName = strcat(path_test_result,'/err_',testCase,'.txt');
        fileID = fopen(fileName,'w');
        fprintf(fileID,'%7s %12s %6s %12s %6s\n','h','L2','CR','ENorm','CR');
        fprintf(fileID,'%6.5f & %12.8f & %6.2f & %12.8f & %6.2f \n',errFileMat);
        fclose(fileID);
    end
else % just for plotting, don't wanna find CR
    
    if savePlot
        
        % w
        nf=nf+1;
        g=figure(nf);
        set(g, 'Visible', 'off');
        plotNXFEM(msh,pa,phi,iPs,nf,whSTD,'withMesh',withMesh,...
                            'title','wh','show',true); % wh
        fileName = strcat(path_test_result,'/',testCase,'_wh_','.png');
        g.PaperUnits = 'inches';
        g.PaperPosition = [0 0 8 6];
        print(fileName,'-dpng','-r0');
        close(g);
        
        % u
        nf=nf+1;
        g=figure(nf);
        set(g, 'Visible', 'off');
        plotNXFEM(msh,pa,phi,iPs,nf,uhSTD,'withMesh',withMesh,...
                            'title','uh','show',true); % wh
        fileName = strcat(path_test_result,'/',testCase,'_uh_','.png');
        g.PaperUnits = 'inches';
        g.PaperPosition = [0 0 8 6];
        print(fileName,'-dpng','-r0');
        close(g);
        
        % v
        nf=nf+1;
        g=figure(nf);
        set(g, 'Visible', 'off');
        plotNXFEM(msh,pa,phi,iPs,nf,vhSTD,'withMesh',withMesh,...
                            'title','vh','show',true); % wh
        fileName = strcat(path_test_result,'/',testCase,'_vh_','.png');
        g.PaperUnits = 'inches';
        g.PaperPosition = [0 0 8 6];
        print(fileName,'-dpng','-r0');
        close(g);
    end
    
    if showPlot
        nf = plotNXFEM(msh,pa,phi,iPs,nf,whSTD,'withMesh',false,'title','wh','dim',3); % wh
        nf = plotNXFEM(msh,pa,phi,iPs,nf,wExSTD,'withMesh',false,'title','wex','dim',3); % wex
        nf = plotNXFEM(msh,pa,phi,iPs,nf,uhSTD,'withMesh',false,'title','uh','dim',3); % uh
        nf = plotNXFEM(msh,pa,phi,iPs,nf,uExSTD,'withMesh',false,'title','uex','dim',3); % uex
        nf = plotNXFEM(msh,pa,phi,iPs,nf,vhSTD,'withMesh',false,'title','vh','dim',3); % vh
        nf = plotNXFEM(msh,pa,phi,iPs,nf,vExSTD,'withMesh',false,'title','vex','dim',3); % vex
    end
end

close all; % close all figures in this test








