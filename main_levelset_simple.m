%% =======================================================================
% This file is used to test with a very simple level set function
% ------------------------------------------------------------------------
% PURPOSE: Coding level set + NOT YET couple with NXFEM (only phi)
% MODEL: model_levelset_*
% ------------------------------------------------------------------------
tic


%% add path of functions
addpath(genpath('func')); % add all necessary functions

% nseg_array = [37, 57, 77];
nseg_array = 57;
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
velo = model.velo; % Velocity (grad of potential in other cases)



%% setting
useFFmesh = 0; % use freefem++ mesh or not?
    reguMesh = 0; % use regular mesh or not? (for matlab mesh)
pa.smallCut = 0; % ignore small-support basis (1=ignore,0=no)
    pa.tH = 10; % to find the small support using (20) or (21) in arnold 2008
useFMM = 1; % use fast marching method or not (mshdist)?
    numUse = 0; % count the number of use of FMM
    alp_FMM = 0.05;
useSUPG = 1; % if 1, need to make more settings
    delEps = 1e-3;
    delSD = 0.5;
showPlot = 0; % show plot or not?
savePlot = 1; % wanna save plot or not?
    pathOption = '_BLOUZA';
    testCase = 'Using ff mesh + wSUPG + wtFMM + thesis (no limit num of use FMM)'; % note for info_and_errors.txt
withMesh = false;

%% choose the machine to run
machine = "blouza"; % options: thi, gia, lehoan, blouza, gaia


% only enable showPlot option on thi's machine
if machine ~= "thi"
    showPlot = 0;
end

for iii=1:size(nseg_array,2)
disp("Running...");
%% mesh setting up
nSeg = nseg_array(iii);
% nSeg = iii;
if ~useFFmesh % don't use ff mesh
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

switch machine
    case "thi"
        path_nxfem = '/home/thi/Dropbox/git/nxfem/'; % thi's local machine
        path_phi = strcat(path_nxfem,'mshdist/');
        call_mshdist = strcat({'mshdist'},{' '},{path_phi},'phi'); % run in terminal
    case "gia"
        path_nxfem = '/home/gia/nxfem/'; % gia's local machine
        path_phi = strcat(path_nxfem,'mshdist/');
        call_mshdist = strcat({'mshdist'},{' '},{path_phi},'phi'); % run in terminal
    case "lehoan"
        path_nxfem = '/home/lehoan/git/nxfem/'; % lehoan's local machine
        path_phi = strcat(path_nxfem,'mshdist/');
        call_mshdist = strcat({'mshdist'},{' '},{path_phi},'phi'); % run in terminal
    case "blouza"
        path_nxfem = '/users/home/blouza/thi/nxfem/'; % blouza's local machine
        path_phi = strcat(path_nxfem,'mshdist/');
        call_mshdist = strcat({'/users/home/blouza/MshDist/build/mshdist'},{' '},{path_phi},'phi'); % run in terminal
    case "gaia" % CHECK LATER!!!!
        path_nxfem = '/users/dinh/nxfem/'; % only on gaia machine
        path_phi = strcat(path_nxfem,'mshdist/');
%         call_mshdist = strcat({'mshdist'},{' '},{path_phi},'phi'); % run in terminal
end
call_mshdist = cell2mat(call_mshdist);



%% write to files .mesh and .sol to use toolbox mshdist
if useFMM
    disp("Write to phi.mesh & phi.sol...");
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
    disp("Creating folder to save plots...");
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
    disp("Show plot of phi...");
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
    disp("Saving plots...");
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

veloOld = @(x,y,sub) velo(x,y,0,sub);



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

disp("Starting loop...");
%% loops
for ns = 1:maxStep
    %% update velo
    t = t+dt;
    
    % if u doesn't depend on time
    % veloNew = veloOld;
    
    % if u depends on time
    veloNew = @(x,y,sub) velo(x,y,t,sub);

    
    
    % CFL
    disp("CFL...");
    CFL(ns) = getCFL(msh, veloNew, dt, hTmax);
    
   
    %% boundary condition
    % Arnold Book p.221, due to velo=0 on boundary, we don't need BCs for phi
    
    
    
    %% get del_T
    if useSUPG
        delOld = getDellsT(msh,veloOld,delEps,delSD); % Arnold's book p.223
        delNew = getDellsT(msh,veloNew,delEps,delSD);
    else
        delOld = zeros(1,size(msh.t,2)); % without SUPG
        delNew = delOld;
    end

   


    %% stiffness matrix
    Enew = getMEls(msh,pa,veloNew,delNew,1);
    Hnew = getMHls(msh,pa,veloNew,delNew,dt*0.5);
    
    % if u doesn't depend on time
    % Aphi = Enew + Hnew;
    
    % If u depends on t
    mI = speye(msh.nStd); % identity matrix
    Aphi = mI + Enew^(-1)*Hnew;
    
    
    
    %% load vector
    
    % If u doesn't depend on t
    % AFphi = Eij - Hij;
    
    % If u depends on t
    Eold = getMEls(msh,pa,veloOld,delOld,1);
    Hold = getMHls(msh,pa,veloOld,delOld,dt*0.5);
    AFphi = mI - Eold^(-1)*Hold;
    
    phi = phi'; % row to column
    Fphi = AFphi*phi;
    
    
    
    %% get new phi
    phi = Aphi\Fphi; % update phi
    phi = phi'; % column to row
    
    
    
    %% Update velocity
    disp("Updating velocity...");
    veloOld = veloNew;

    
    
    norm_gradphi = getNormL2Gstd(msh,phi); % ||gradPhi||_L2
    %% mshdist
%     if useFMM && abs(1-norm_gradphi) > alp_FMM && numUse <=1
    if useFMM && abs(1-norm_gradphi) > alp_FMM
        disp("Use FMM...");
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
    disp("Saving error files...");
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

toc

%% get deltaT: 1 x nTs
function delT = getDellsT(msh,vel_t,eps,SD)
    % This function is only to be used type of velocity = (x,y,sub)
    % This function follows formula presented in Arnold's book, page 223 (and also cf. (7.14))
    % If in future, we need it more than once, I will put it in a separated file
    % Input: - msh.hT : h on each triangle : 1 x nTs
    %        - velo at each time step t: (x,y,sub)
    %        - msh.p to get all x,y coordinates
    %        - a control number SD = O(1) (cf. (7.14))
    %        - given small eps > 0 (cf. 7.14), at page 223, he took eps=1e-3
    % Output: a vector 1 x nTs

    x = msh.p(1,:); y = msh.p(2,:);

    % Take velo at each triangle
    veloPu = vel_t(x,y,1); % 1 x nPs
    veloPv = vel_t(x,y,2); % 1 x nPs
    veloTu = veloPu(msh.t(1:3,:)); % 3 x nTs
    veloTumax = max(abs(veloTu),[],1); % 1 x nTs
    veloTv = veloPv(msh.t(1:3,:)); % 3 x nTs
    veloTvmax = max(abs(veloTv),[],1); % 1 x nTs
    normVeloT = max(veloTumax,veloTvmax); % 1 x nTs

    normVeloT = max(eps,normVeloT); % 1 x nTs

    delT = SD*msh.hT./normVeloT;
end

function val = getCFL(msh, velo, dt, dx)
    x = msh.p(1,:); y = msh.p(2,:);
    
    % Take velo at each triangle
    veloPu = velo(x,y,1); % 1 x nPs
    veloPv = velo(x,y,2); % 1 x nPs
    veloTu = veloPu(msh.t(1:3,:)); % 3 x nTs
    veloTumax = max(abs(veloTu),[],1); % 1 x nTs
    veloTv = veloPv(msh.t(1:3,:)); % 3 x nTs
    veloTvmax = max(abs(veloTv),[],1); % 1 x nTs
    normVeloT = max(veloTumax,veloTvmax); % 1 x nTs
    
    val = max(normVeloT)*dt/dx;
end