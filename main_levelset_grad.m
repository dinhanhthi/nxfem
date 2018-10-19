%% =======================================================================
% This file is used to test with a very simple level set function
% The same with main_levelset_simple BUT ONLY FOR GRAD
% grad(u) in this test case = u in main_levelset_simple
% ------------------------------------------------------------------------
% PURPOSE: Coding level set + NOT YET couple with NXFEM (only phi)
% MODEL: model_levelset_*
% ------------------------------------------------------------------------
% RESULT: Cannot verify because vh in Vh, we cannot choose a v by a normal
% expression to obtain grad of v as usual (grad of v in Vh is totally
% different!!!)


%% add path of functions
addpath(genpath('func')); % add all necessary functions

% nseg_array = [37, 57, 77];
nseg_array = 27;
% nseg_array = zeros(1,3); % get from ffpp 



%% Fixed parameters
pa.degP1D = 3; % Gaussian quadrature points in 1D (polinomial functions)
pa.degP2D = 4; % Gaussian quadrature points in 2D (polinomial functions)
pa.degN = 8; % Gaussian quadrature points in 2D (non-polynomial functions)
% degree-#OfPoints : 1-1, 2quick -3, 3-4, 4-6, 5-7, 6-12, 7-13, 8-16, 9-19, 10-25, 11-27, 12-33
pa.tol = eps(1e3); % tolerance, 1e-14


%% models
model = model_levelset_x; % Becker's test case with interface: x=x_0
% model = model_levelset_vortex;  % Niklas' test case
GeoDom = model.domain(); % domain
veloGrad = model.veloGrad; % Velocity (grad of potential in other cases)



%% setting
useFFmesh = 0; % use freefem++ mesh or not?
    reguMesh = 0; % use regular mesh or not? (for matlab mesh)
pa.smallCut = 0; % ignore small-support basis (1=ignore,0=no)
    pa.tH = 10; % to find the small support using (20) or (21) in arnold 2008
useFMM = 0; % use fast marching method or not (mshdist)?
    numUse = 0; % count the number of use of FMM
    alp_FMM = 0.1;
useSUPG = 0; %nseg_array = 17; if 1, need to make more settings
    delEps = 1e-3;
    delSD = 0.5;
showPlot = 1; % show plot or not?
savePlot = 0; % wanna save plot or not?
    pathOption = '_FF_new';
    testCase = 'Using ff mesh + wSUPG + wtFMM + thesis (no limit num of use FMM)'; % note for info_and_errors.txt
withMesh = false;


for iii=1:size(nseg_array,2)
disp("Running...");
%% mesh setting up
nSeg = nseg_array(iii);
% nSeg = iii;
if ~useFFmesh
    disp("Use matlab mesh...");
    if ~reguMesh % not regular mesh?
        hEdgeMax = 2/nSeg;
        [points,edges,triangles] = initmesh(GeoDom,'hmax',hEdgeMax); %irregular
    else
        [points,edges,triangles] = poimesh(GeoDom,nSeg,nSeg); % regular
    end
else
    disp("Use FreeFem mesh...");
    if ~reguMesh
        file = strcat('Th_irregular_new',num2str(iii-1),'.msh');
    else
        file = strcat('Th_regular',num2str(iii-1),'.msh');
    end
    [points,edges,triangles] = getMeshFromFF(file);
end
msh.p = points; msh.t = triangles; msh.e = edges; % save to msh
x = points(1,:); % x-coordinate of points
y = points(2,:); % y-coordinate of points
% diameter (longest side) of each triangle: 1 x nTs
msh.hT = getDiam(msh); % 1 x number of triangles
msh.hTmax = max(msh.hT); % maximum of all diameters
hTmax = msh.hTmax;



%% Level set function (INITIAL)
phi = model.defPhi(x,y,pa); % 1 x number of points (row array)
% phi(abs(phi)<pa.tol) = 0; % find phi which are very small (~0) and set to 0



%% create command to run mshdist outside matlab
path_nxfem = '/home/thi/Dropbox/git/nxfem/'; % thi's local machine
% path_nxfem = '/users/dinh/nxfem/'; % only on gaia machine
path_phi = strcat(path_nxfem,'mshdist/');
call_mshdist = strcat({'mshdist'},{' '},{path_phi},'phi'); % run in terminal
call_mshdist = cell2mat(call_mshdist);



%% write to files .mesh and .sol to use toolbox mshdist
if useFMM
    mshdist_w_mesh(msh,path_phi,'phi'); % export to phi.mesh
    mshdist_w_sol(msh,phi,path_phi,'phi'); % export to phi.sol
end



%% setting up for the evolution problem
pA = model.pa;
pA = pA();
t=0;
dt = pA.dt;
Tmax = pA.Tmax;
maxStep = Tmax/dt;



%% create folder of test files results
if savePlot
    if reguMesh
       path_regu = 'regu_';
    else
        path_regu = 'irregu_';
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
    path_test_result = strcat(path_nxfem,'results/level_set/',...
                                path_regu,num2str(nSeg),'_',num2str(dt),...
                                    path_useSUPG,path_useFMM,pathOption);
    path_test_remove = strcat({'rm -r'},{' '},{path_test_result}); % in case of duplicated folder
    path_test_remove = cell2mat(path_test_remove);
    system(path_test_remove);
    path_test_create = strcat({'mkdir'},{' '},{path_test_result}); % crfeate a new folder
    path_test_create = cell2mat(path_test_create);
    system(path_test_create);
end



%% initial phi and plot
tris = getTriangles(phi,msh,pa); % tris has 3 factors (structure var)
CTs=tris.CTs;
msh.nStd = size(points,2); % number of standard nodes
CT = getInfoCTs(CTs,phi,msh,pa); % CT has many factors (structure var)
iPs = CT.iPs;

nf = 0; % reset every loop to be sure uh, vh plotted on the same figure
titlePlot = strcat('t= ',num2str(t));

if showPlot
    disp("Plotting phi...");
    plotNXFEM(msh,iPs,nf,phi,'withMesh',withMesh,'title',titlePlot,'dim',2,'export',false); % phi
    nCTs = size(iPs,3);
    hold on
    for it=1:nCTs
        plot(iPs(1,:,it),iPs(2,:,it),'-r','LineWidth',1);
        hold on
    end
    hold off
end

% if useFMM
%     mshdist_w_sol(msh,phi,path_phi,'phi'); % export to phi.sol
%     system(call_mshdist); % run 'mshdist file/to/phi' (redistancing)
%     phi = mshdist_r_sol(phi,path_phi,'phi'); % update phi
% end

if savePlot
    f=figure('visible','off');
    pdeplot(points,edges,triangles,'XYData',phi,'Title',titlePlot);
    hold on
    if withMesh
        pdemesh(points,edges,triangles);
    end
    nCTs = size(CTs,2);
    for tri=1:nCTs
        plot(iPs(1,:,tri),iPs(2,:,tri),'-r','LineWidth',1);
        hold on
    end
    hold off
    fileName = strcat(path_test_result,'/ls_',num2str(0),'_','t_',num2str(t),'.png');
    
    % change size of images
    oldpaperunits = get(gcf,'PaperUnits');
    oldpaperpos = get(gcf,'PaperPosition');
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 6];
    print(fileName,'-dpng','-r0');
    
%     print(fileName,'-dpng'); % old
    close(f);
end



%% errors (cf. page 224, Arnold book)
phi0h = phi;
if ~isempty(CTs)
    normPhih0oGh0 = getNormL2oGstd(msh,pa,CTs,CT,phi0h); % initial phi
    normPhionGh0 = getNormL2foGh(msh,pa,CTs,iPs,model.defPhi); % ||phi^0||_L2(Gam_h^0)
%     display(normPhionGh0);
end

vhOld = zeros(msh.nStd, 1); % initial velocity



% %% Plot velo with mesh
% f=figure('visible','off');
% velo_x = veloOld(x,y,1);
% velo_y = veloOld(x,y,2);
% set(gca,'FontSize',24) % Creates an axes and sets its FontSize to 18
% quiver(x,y,velo_x,velo_y);
% set(gca,'FontSize',24);
% fileName = strcat('results/level_set/velo_t0_nSeg_17.eps');
% oldpaperunits = get(gcf,'PaperUnits');
% oldpaperpos = get(gcf,'PaperPosition');
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 8 6];
% print(fileName,'-depsc','-r0');
% close(f);


CFL = zeros(1,maxStep); % store CFL value

disp("Starting the loop...");
%% loops
for ns = 1:maxStep
    disp("-----------------------------");
    Xdisp = ['step = ', num2str(ns)];
    disp(Xdisp);
    
    %% update velo
    t = t+dt;
    
    % if u doesn't depend on time
    % veloNew = veloOld;
    
    % if u depends on time
%     vhSTD = veloGrad(x,y,t);
    vhSTD = vhOld;

    
    
    % CFL
%     CFL(ns) = getCFL(msh, vhSTD, dt, hTmax);
    
   
    %% boundary condition
    % Arnold Book p.221, due to velo=0 on boundary, we don't need BCs for phi
    
    
    
    %% get del_T
    if useSUPG
        delOld = getDellsT(msh,vhOld,delEps,delSD); % Arnold's book p.223
        delNew = getDellsT(msh,vhSTD,delEps,delSD);
    else
        delOld = zeros(1,size(msh.t,2)); % without SUPG
        delNew = delOld;
    end

   


    %% stiffness matrix
    Enew = getMEls_gP(msh,pa,vhSTD,delNew,1);
    Hnew = getMHls_gP(msh,pa,vhSTD,delNew,dt*0.5);
    mI = speye(msh.nStd); % identity matrix
    Aphi = mI + Enew^(-1)*Hnew;
    
    
    
    %% load vector
    Eold = getMEls_gP(msh,pa,vhOld,delOld,1);
    Hold = getMHls_gP(msh,pa,vhOld,delOld,dt*0.5);
    AFphi = mI - Eold^(-1)*Hold;
    phi = phi'; % row to column
    Fphi = AFphi*phi;
    
    
    
    %% get new phi
    phi = Aphi\Fphi; % update phi
    phi = phi'; % column to row
    
    
    
    %% Update velocity
    disp("Updating v...");
    vhOld = vhSTD;

    
    
    norm_gradphi = getNormL2Gstd(msh,phi); % ||gradPhi||_L2
    %% mshdist
%     if useFMM && abs(1-norm_gradphi) > alp_FMM && numUse <=1
    if useFMM && abs(1-norm_gradphi) > alp_FMM
        disp("FMM...");
        mshdist_w_sol(msh,phi,path_phi,'phi'); % export to phi.sol
        system(call_mshdist); % run 'mshdist file/to/phi' (redistancing)
        phi = mshdist_r_sol(phi,path_phi,'phi'); % update phi
        numUse = numUse + 1;
    end
    
    
    
    
    %% triangle and mesh info (of new phi)
    % Just for plotting and saving the plots
    tris = getTriangles(phi,msh,pa); % tris has 3 factors (structure var)
    CTs = tris.CTs;
    CT = getInfoCTs(CTs,phi,msh,pa); % CT has many factors (structure var)
    iPs = CT.iPs;
    
    
    
    
    %% plot phi
%     abc = waitforbuttonpress; % wait for click
    disp("Plotting phi");
    titlePlot = strcat('t= ',num2str(t));
    if showPlot
        plotNXFEM(msh,iPs,nf,phi,'withMesh',withMesh,'title',titlePlot,'dim',2,'export',false); % phi
        nCTs = size(iPs,3);
        hold on
        for it=1:nCTs
            plot(iPs(1,:,it),iPs(2,:,it),'-r','LineWidth',1);
            hold on
        end
        hold off
    end
    
    
    
%     %% ||grad phi||
%     norm_gradphi = getNormL2Gstd(msh,phi)
    
    
    
    %% save plot to file
    if savePlot
        disp("Saving plot...");
        if ((mod(ns,4)==0)&&(ns < maxStep-5)) || (ns >= maxStep-5)
            f=figure('visible','off');
            pdeplot(points,edges,triangles,'XYData',phi,'Title',titlePlot);
            hold on
            if withMesh
                pdemesh(points,edges,triangles);
            end
            nCTs = size(CTs,2);
            for tri=1:nCTs
                plot(iPs(1,:,tri),iPs(2,:,tri),'-r','LineWidth',1);
                hold on
            end
            hold off
            fileName = strcat(path_test_result,'/ls_',num2str(ns),'_','t_',num2str(t),'.png');

            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 8 6];
            
            % change size of images
            print(fileName,'-dpng','-r0');
    
%             print(fileName,'-dpng'); % old
            close(f); 
        end
    end
    
end % for ns



%% errors
normPhihNoGhN = getNormL2oGstd(msh,pa,CTs,CT,phi); % ||phi_h^N||_{L2(Gam_h^N)}
err = phi - phi0h;
normErrPhihN = getNormL2std(msh,pa,err); % ||phi_h^N-phi_h^0||_{L2(Omg)}
normPhionGhN = getNormL2foGh(msh,pa,CTs,iPs,model.defPhi); % ||phi^0||_L2(Gam_h^N)
normPhih0PhionOmg = getNormL2fhf(msh,pa,phi0h,model.defPhi); % ||phi_h^0 - phi^0||_L2(Omg)
normPhihNPhionOmg = getNormL2fhf(msh,pa,phi,model.defPhi); % ||phi_h^N - phi^0||_L2(Omg)



%% display results
% fprintf('hTmax: %0.7f,\n',msh.hTmax);
% % fprintf('||phi_h^0||_{L2(Gam_h^0)}: %0.10f,\n',normPhih0oGh0);
% % fprintf('||phi_h^N||_{L2(Gam_h^N)}: %0.10f,\n',normPhihNoGhN);
% fprintf('||phi_h^N-phi_h^0||_{L2(Omg)}: %0.10f\n',normErrPhihN);
% fprintf('||phi^0||_L2(Gam_h^N): %0.10f\n',normPhionGhN);
% fprintf('||phi^0||_L2(Gam_h^0): %0.10f\n',normPhionGh0);
% fprintf('||phi_h^0-phi^0||_{L2(Omg)}: %0.10f\n',normPhih0PhionOmg);
% fprintf('||phi_h^N-phi^0||_{L2(Omg)}: %0.10f\n',normPhihNPhionOmg);



%% save error result to file
if savePlot
    disp("Saving errors...");
    fileName = strcat(path_test_result,'/info_and_errors_',num2str(nSeg),'.txt');
    fileID = fopen(fileName,'w');
        fprintf(fileID,'%s,\n',testCase);
        fprintf(fileID,'Model: %s,\n',model.name);
        fprintf(fileID,'nSeg: %d,\n',nSeg);
        fprintf(fileID,'hTmax: %0.10f,\n',hTmax);
        fprintf(fileID,'dt: %f,\n',dt);
        fprintf(fileID,'Regular mesh: %d,\n',reguMesh);
        fprintf(fileID,'Use small-cut: %d,\n',pa.smallCut);
        fprintf(fileID,'Use FMM: %d,\n',useFMM);
        fprintf(fileID,'al_FMM: %f,\n',alp_FMM);
        fprintf(fileID,'numUseFMM: %d,\n',numUse);
        fprintf(fileID,'useSUPG: %d,\n',useSUPG);
        fprintf(fileID,'delEps: %f,\n',delEps);
        fprintf(fileID,'delSD: %f,\n',delSD);
        fprintf(fileID,'CFL: ');
        for iclf = 1:size(CFL,2)
            fprintf(fileID,'%0.2f, ',CFL(1,iclf));
        end
        fprintf(fileID,'\n');
        fprintf(fileID,'\n||phi_h^0||_{L2(Gam_h^0)}: %0.10f,\n',normPhih0oGh0);
        fprintf(fileID,'||phi_h^N||_{L2(Gam_h^N)}: %0.10f,\n',normPhihNoGhN);
        fprintf(fileID,'||phi_h^N-phi_h^0||_{L2(Omg)}: %0.10f.\n',normErrPhihN);
        fprintf(fileID,'||phi^0||_L2(Gam_h^N): %0.10f\n',normPhionGhN);
        fprintf(fileID,'||phi^0||_L2(Gam_h^0): %0.10f\n',normPhionGh0);
        fprintf(fileID,'||phi_h^0-phi^0||_{L2(Omg)}: %0.10f.\n',normPhih0PhionOmg);
        fprintf(fileID,'||phi_h^N-phi^0||_{L2(Omg)}: %0.10f.\n',normPhihNPhionOmg);
    fclose(fileID);
end
    
end

%% get deltaT: 1 x nTs
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

% function val = getCFL(msh, velo, dt, dx)
%     x = msh.p(1,:); y = msh.p(2,:);
%     
%     % Take velo at each triangle
%     veloPu = velo(x,y,1); % 1 x nPs
%     veloPv = velo(x,y,2); % 1 x nPs
%     veloTu = veloPu(msh.t(1:3,:)); % 3 x nTs
%     veloTumax = max(abs(veloTu),[],1); % 1 x nTs
%     veloTv = veloPv(msh.t(1:3,:)); % 3 x nTs
%     veloTvmax = max(abs(veloTv),[],1); % 1 x nTs
%     normVeloT = max(veloTumax,veloTvmax); % 1 x nTs
%     
%     val = max(normVeloT)*dt/dx;
% end