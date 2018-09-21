%% ========================================================
% INTRODUCTION
% =========================================================
% NXFEM based on Hansbo's paper.
% This file is used for both Sinha's models and Barrau's model
% =========================================================



%% add path of functions
addpath(genpath('func')); % add all necessary functions



%% Fixed parameters
pa.degP1D = 3; % Gaussian quadrature points in 1D (for polinomial functions)
pa.degP2D = 4; % Gaussian quadrature points in 2D (for polinomial functions)
pa.degN = 8; % Gaussian quadrature points in 2D (for non-polynomial functions)
% degree-#OfPoints : 1-1, 2-3, 3-4, 4-6, 5-7, 6-12, 7-13, 8-16, 9-19, 10-25, 11-27, 12-33
pa.tol = eps(1e3); % tolerance, 1e-14



%% model 
mdl = 1; 
switch mdl
    case 1 % sinha (cf. model_sinha.m)
        model=model_sinha; % file sinha.m
        pa.kk1 = 1; pa.kk2 = 0.5; % diffusion coefficients (kk1 must be = 2kk2)
        pa.xi = 1.0; % DON'T CHANGE!!!
    case 2 % barrau (p.44) STRAIGHT LINE [-1,1]x[-1,1]
        model=model_barrau_x; % file barrau_x.m
        pa.kk1 = 1; pa.kk2 = 1000; % for barrau case
        pa.xi = 0.91; % only used for barrau's model, between (-1,1)
    case 3 % barrau_r (p.47) PART OF CIRCLE [0,1]x[0,1]
        model=model_barrau_r; % file barrau_r.m
        pa.kk1 = 1; pa.kk2 = 1000; % for barrau case
%         pa.kk1 = 1; pa.kk2=1; % testing
        pa.xi = 0.75; % only used for barrau's model, between (0,1)
    case 4 % barrau_c (p.113) CIRCLE INSIDE DOMAIN [-1,1]x[-1,1]
        model=model_barrau_c;
        pa.kk1 = 1; pa.kk2 = 1000;
        pa.xi = 0.71;
    case 5 % only w in article1
        model=model_w;
%         pa.kk1=1; pa.kk2=100;
        pa.kk1=1; pa.kk2=1;
        pa.xi=0.6;
        pa.lamSys=1;
    case 6 % test case in Boiveau's thesis (cf. 5.5)
        % note that, this is for nonsysmmetric + penalty free NXFEM method
        model = model_boiveau;
        pa.kk1=1; pa.kk2 = 100;        
%         pa.kk1=10; pa.kk2 = 1;
        pa.xi=0.3;
end
GeoDom = model.domain(); % domain



%% setting
findCR = 1; % wanna find convergence rate? 1 or 0
useFilesNotPlot = 0; % 1=export figures to files (and without plotting them)
pa.smallCut = 1; % ignore small-support basis (1=ignore,0=no)
pa.tH = 100; % to find the small support using (20) or (21) in arnold 2008
pa.lamH = 1e4; % penalty coefficient



%% ghost penalty
pa.useGP = 1; % wanna use ghost penalty term?
pa.gam1 = 1e-9; % parameter for 1st term
pa.gam2 = 1e-9; % parameter for 2nd term



%% more settings
if findCR==1 % wanna  find the convergence rate
    nStep = 4; % number of interations
    pa.reguMesh = 0; % use regular mesh or not?
    wannaPlot = 0; % DON'T CHANGE!
    useDOF = 0; % use dofs to find convergence rate? (default = h)
else % should be used only for plotting
    nStep = 1; % shouldn't change
    pa.reguMesh = 0; % use regular mesh or not?
    nSeg = 51; % ONLY FOR nStep=1;
    wannaPlot = 1; % wanna plot or not? (JUST FOR nStep=1)
    % Go to the end of this file to change what you wanna plot
end



%% ========================================================
% VARIABLES
% =========================================================
hEdgeMax = zeros(nStep,1); % maximum edge size (edge of domain)
%% structure variables
% mesh msh, errors err, solutions sol, plotting components pplot


%% ========================================================
% MAIN PROCESS
% =========================================================
if nStep>1
    numSeg = [41 51 71 91];
    % numSeg = [81 91 101 121];
else
    numSeg(1)=nSeg; 
end
for z=nStep:-1:1
%     if nStep>1
% %         numSeg(z) = 2^(z-1)*10+1; % 11,21,41,81
%         numSeg(z) = 2^z*10+1; % 21,41,81,161
%     else % nStep=1
%         numSeg(z) = nSeg;
%     end
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
    tx=zeros(1,nStep); ty=zeros(2,nStep);
    for i=1:nStep
       tx(i)=msh(i).hTmax;
       ty(1,i)=err(i).L2;
       ty(2,i)=err(i).ENorm;
    end
    tmp = polyfit(log(tx),log(ty(1,:)),1);
    order.L2 = tmp(1);
    tmp = polyfit(log(tx),log(ty(2,:)),1);
    order.ENorm = tmp(1);
    
    
    %% Figure of convergence rate
    % Only for error in L2 and ENorm
    % See the file results\errPlot.jpg
    if useFilesNotPlot
        f=figure('visible','off');
        plot(log(tx),log(ty(1,:)),'-r.',...
            log(tx),log(ty(2,:)),'-b.');
        legend('L2','ENorm');
        xlabel('log(h)'); 
        ylabel('log(error)');
        print -djpeg results\main_err.jpg;
        close(f);
        
        %% EXPORT errors TO FILES
        % see file results\errors.txt
        
        cr = zeros(2,nStep);
        for i=2:nStep
           cr(1,i) = log(ty(1,i)/ty(1,i-1))/log(tx(i)/tx(i-1)); % L2
           cr(2,i) = log(ty(2,i)/ty(2,i-1))/log(tx(i)/tx(i-1)); % ENorm
        end
        errFileMat=[tx;ty(1,:);cr(1,:);ty(2,:);cr(2,:)];
        fileID = fopen('results\main_err.txt','w');
        fprintf(fileID,'Model: %2.0f \n',mdl);
        fprintf(fileID,'%7s %12s %6s %12s %6s\n','h','L2','CR','ENorm','CR');
        fprintf(fileID,'%6.5f & %12.8f & %6.2f & %12.8f & %6.2f \n',errFileMat);
        fclose(fileID);
    else
        plot(log(tx),log(ty(1,:)),'-r.',...
            log(tx),log(ty(2,:)),'-b.');
        legend('L2','ENorm');
        xlabel('log(h)'); 
        ylabel('log(error)');
    end
    
    
    % display results
    fprintf('regular: %d,\t',pa.reguMesh);
    fprintf('numSeg: %0.0f,\t',numSeg(nStep));
    fprintf('DOFs: %0.0f\n',msh(nStep).ndof);
    fprintf('L2 order: %0.15f\n',order.L2);
    fprintf('ENorm order: %0.15f\n',order.ENorm);
else % only 1 step
    fprintf('numSeg: %0.0f,\t',numSeg(nStep));
    fprintf('DOFs: %0.0f\n',msh(nStep).ndof);
    fprintf('L2 err: %0.18f\n',err(nStep).L2);
    fprintf('ENorm err: %0.18f\n',err(nStep).ENorm);
end


%% ========================================================
% PLOTTING
% =========================================================
if wannaPlot==1
    nf = 0;
    % nf = plotNXFEM(msh,pplot.iPs,nf,'eleLabel','off','nodeLabel','on'); % only mesh
%     nf = plotNXFEM(msh,pplot.iPs,nf,sol(1).Vh,'withMesh',false,'title','uh','dim',2,'withGamh',true); % uh
%     nf = plotNXFEM(msh,pplot.iPs,nf,sol(1).ex,'withMesh',false,'title','uex','dim',2,'withGamh',true); % uex
    nf = plotNXFEM(msh,pplot.iPs,nf,sol(1).ex,'withMesh',false,'title','uex','dim',3,'export',false); % uex
    nf = plotNXFEM(msh,pplot.iPs,nf,sol(1).Vh,'withMesh',false,'title','uh','dim',3); % uh
    % nf = plotNXFEM(msh,pplot.iPs,nf,phi{1},'withMesh',false,'title','level set','dim',3,'export',false); % level set function
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
