%% ========================================================
% INTRODUCTION
% =========================================================
% NXFEM based on Hansbo's paper.
% This file is used for both Sinha's models and Barrau's model
% =========================================================


%% ========================================================
% FIXED PARAMERTERS
% shouldn't change
% =========================================================
% pa = all parameters
pa.degP1D = 3; % Gaussian quadrature points in 1D (for polinomial functions)
pa.degP2D = 4; % Gaussian quadrature points in 2D (for polinomial functions)
pa.degN = 8; % Gaussian quadrature points in 2D (for non-polynomial functions)
% degree-#OfPoints : 1-1, 2-3, 3-4, 4-6, 5-7, 6-12,
%                    7-13, 8-16, 9-19, 10-25, 11-27, 12-33
pa.tol = eps(1e3); % tolerance, 1e-14


%% ========================================================
% INPUT REGION
% =========================================================
% available models: sinha, barrau (2 cases)
mdl = 2; % 1 for sinha, 2 for barrau with x, 3 for barrau with circle
findCR = 1; % wanna find convergence rate? 1 or 0
useFilesNotPlot = 0; % 1=export figures to files (and without plotting them)
pa.smallCut = 0; % ignore small-support basis (1=ignore,0=no)
pa.tH = 10; % to find the small support using (20) or (21) in arnold 2008

%------------------------------------
% PARAMETERS
%------------------------------------
pa.lambdaH = 1e8; % penalty coefficient
pa.kap1per2 = 0; % kap1=kap2=1/2 or not? "1"=yes, "0"=no

%------------------------------------
%% GHOST PENALTY
%------------------------------------
pa.useGP = 1; % wanna use ghost penalty term?
% There are some bugs (fix later)
pa.gam1 = 1e-1; % parameter for 1st term
pa.gam2 = 1e-1 ; % parameter for 2nd term

%------------------------------------
% SETTINGS
%------------------------------------
if findCR==1 % wanna  find the convergence rate
    nStep = 4; % number of interations
    pa.reguMesh = 0; % use regular mesh or not?
    wannaPlot = 0; % DON'T CHANGE!
    useDOF = 0; % use dofs to find convergence rate? (default = h)
else % should be used only for plotting
    nStep = 1; % shouldn't change
    pa.reguMesh = 0; % use regular mesh or not?
    nSeg = 25; % ONLY FOR nStep=1;
    wannaPlot = 1; % wanna plot or not? (JUST FOR nStep=1)
    % Go to the end of this file to change what you wanna plot
end

%------------------------------------
%% MODELS
%------------------------------------
switch mdl
    case 1 % sinha
        model=sinha; % file sinha.m
        pa.kk1 = 1; pa.kk2 = 0.5; % diffusion coefficients (kk1 must be = 2kk2)
        pa.xi = 1.0; % DON'T CHANGE!!!
    case 2 % barrau (p.44) STRAIGHT LINE [-1,1]x[-1,1]
        model=barrau_x; % file barrau_x.m
        pa.kk1 = 1; pa.kk2 = 1000; % for barrau case
%         pa.kk1 = 1; pa.kk2 = 100;
%         pa.xi = 0.91; % only used for barrau's model, between (-1,1)
        pa.xi=0.3; % just for testing
    case 3 % barrau_r (p.47) PART OF CIRCLE [0,1]x[0,1]
        model=barrau_r; % file barrau_r.m
%         pa.kk1 = 1; pa.kk2 = 1000; % for barrau case
        pa.kk1 = 1; pa.kk2 = 10; % for barrau case
        pa.xi = 0.75; % only used for barrau's model, between (0,1)
    case 4 % barrau_c (p.113) CIRCLE INSIDE DOMAIN [-1,1]x[-1,1]
        model=barrau_c;
%         pa.kk1 = 1; pa.kk2 = 1000;
        pa.kk1=2; pa.kk2=20;
        pa.xi = 0.75;
end
GeoDom = model.domain(); % domain


%% ========================================================
% VARIABLES
% =========================================================
numSeg = zeros(nStep,1); % number of segments
hEdgeMax = zeros(nStep,1); % maximum edge size (edge of domain)
%% structure variables
% mesh msh, errors err, solutions sol, plotting coponents pplot


%% ========================================================
% MAIN PROCESS
% =========================================================
% for z=1:nStep
for z=nStep:-1:1
    if nStep>1
%         numSeg(z) = 2^z*10+1; % 21,41,81,161
        numSeg(z) = 2^z*10+1; % 11,21,41,81
%         numSeg(z) = 20*(z+2)+1; % 61,81,101,121 
%         numSeg = [41,81,101,121];
    else % nStep=1
        numSeg(z) = nSeg;
    end
    hEdgeMax(z) = 2/numSeg(z); % maximum edge size (edge of domain)
    [msh(z),err(z),sol(z),pplot(z),phi{z}] =...
                        main_eachStep(hEdgeMax(z),GeoDom,pa,model);
end


%% ========================================================
% CONVERGENCE RATE
% display convergence rate
% =========================================================
if nStep>1 % finding convergence rate
    
    %% finding the convergence rate
    tx=zeros(1,nStep); ty=zeros(4,nStep);
    for i=1:nStep
       tx(i)=msh(i).hTmax;
       ty(1,i)=err(i).L2;
%        ty(2,i)=err(i).H1;
%        ty(3,i)=err(i).L2G;
       ty(4,i)=err(i).ENorm;
    end
    tmp = polyfit(log(tx),log(ty(1,:)),1);
    order.L2 = tmp(1);
%     tmp = polyfit(log(tx),log(ty(2,:)),1);
%     order.H1 = tmp(1);
%     tmp = polyfit(log(tx),log(ty(3,:)),1);
%     order.L2G = tmp(1);
    tmp = polyfit(log(tx),log(ty(4,:)),1);
    order.ENorm = tmp(1);
    
    
    %% Figure of convergence rate
    % Only for error in L2 and ENorm
    % See the file results\errPlot.jpg
    if useFilesNotPlot
        f=figure('visible','off');
        plot(log(tx),log(ty(1,:)),'-r.',...
            log(tx),log(ty(4,:)),'-b.');
        legend('L2','ENorm');
        xlabel('log(h)'); 
        ylabel('log(error)');
        print -djpeg results\errPlot.jpg;
        close(f);
        
        %% EXPORT errors TO FILES
        % see file results\errors.txt
        errFileMat=[tx;ty];
        fileID = fopen('results\errors.txt','w');
        fprintf(fileID,'%0s %9s %9s %7s %11s\n',...
                        'h','L2','L2G','H1','ENorm');
        fprintf(fileID,'%0.6f %2.6f %2.6f %2.6f %2.6f\n',errFileMat);
        fclose(fileID);
    else
        plot(log(tx),log(ty(1,:)),'-r.',...
            log(tx),log(ty(4,:)),'-b.');
        legend('L2','ENorm');
        xlabel('log(h)'); 
        ylabel('log(error)');
    end
    
    
    % display results
    fprintf('regular: %d,\t',pa.reguMesh);
    fprintf('kap1per2: %d,\t',pa.kap1per2);
    fprintf('numSeg: %0.0f,\t',numSeg(nStep));
    fprintf('DOFs: %0.0f\n',msh(nStep).ndof);
    fprintf('L2 order: %0.5f\n',order.L2);
%     fprintf('L2Grad order: %0.5f\n',order.L2G);
%     fprintf('H1 order: %0.5f\n',order.H1);
    fprintf('ENorm order: %0.5f\n',order.ENorm);
else % only 1 step
    fprintf('numSeg: %0.0f,\t',numSeg(nStep));
    fprintf('DOFs: %0.0f\n',msh(nStep).ndof);
    fprintf('L2 err: %0.8f\n',err(nStep).L2);
%     fprintf('L2G err: %0.8f\n',err(nStep).L2G);
%     fprintf('H1 err: %0.8f\n',err(nStep).H1);
    fprintf('ENorm err: %0.8f\n',err(nStep).ENorm);
    fprintf('NJU: %0.5f\n',err(nStep).NJU);
    fprintf('NJGU: %0.5f\n',err(nStep).NJGU);
%     fprintf('etaK: %0.5f\n',err(nStep).etaK);
%     fprintf('etaS: %0.5f\n',err(nStep).etaS);
%     fprintf('zetaS: %0.5f\n',err(nStep).zetaS);
end


%% ========================================================
% PLOTTING
% =========================================================
if wannaPlot==1
    D2=2; D3=3; D2D3=1; % which dimension? (D1D2 for both cases)
    wM=1; wtM=0; % plot with mesh or without mesh?
    On = 'on'; Off = 'off';
    % MESH --------------------------------------
    mshs{1} = 0; % wanna plot mesh? (1 or 0)
    mshs{2} = On; % node labels? (On or Off)
    mshs{3} = On; % element label? (On or Off)
    % LEVEL SET FUNCTION ------------------------
    lsf{2} = 0; % wanna plot level set function? (1 or 0)
    lsf{1} = phi{1}; % level set function
    lsf{3} = D3; % dimension? (D2, D3 or D1)
    lsf{4} = wM; % plot with mesh? (wM or wtM)
    % EXACT SOLUTION ----------------------------
    eS{2} = 1; % wanna plot exact solution? (1 or 0)
    eS{1} = sol(1).ex; % exact solution
    eS{3} = D3; % dimension? (D2, D3 or D1)
    eS{4} = wM; % plot with mesh? (wM or wtM)
    % NUMERICAL SOLUTION ------------------------
    nS{2} = 1; % wanna plot numerical solution? (1 or 0)
    nS{1} = sol(1).Vh; % numerical solution (in standard FEM)
%     nS{1} = sol(1).errSolVh; % error u-uex (in standard FEM)
    nS{3} = D3; % dimension? (D2, D3 or D1)
    nS{4} = wM; % plot with mesh? (wM or wtM)
%     % interface Gam_h
%     gamh = pplot.gamh;
    % DO THE PLOT -------------------------------------
    if useFilesNotPlot % just export to the files
        plotNXFEMex(nS,eS,mshs,lsf,msh);
    else % cannot export to files
        plotNXFEM(nS,eS,mshs,lsf,msh,pplot.iPs);
    end
end

plotGrad = 0; % plot gradient
% errVhEx = sol(1).ex-sol(1).Vh;
if plotGrad
   figure('NumberTitle','on','Name','ex');
   plotGradient(sol(1).ex,msh.t,msh,1,4);
   figure('NumberTitle','off','Name','num')
   plotGradient(sol(1).Vh,msh.t,msh,1,5);
   plotGradient(errVhEx,msh.t,msh,1,5);
end

% plot unit normal vector
% if ~isempty(pplot.iPs)
%     plotUNVwM(pplot.iPs,pplot.uNVCTs,msh);
% end
