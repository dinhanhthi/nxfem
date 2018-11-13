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



%% ========================================================
% INPUT REGION
% =========================================================
% available models: sinha, barrau
mdl = 3; % 1 for sinha, 2 for barrau with x, 3 for barrau with circle
findCR = 1; % wanna find convergence rate? 1 or 0
% Parameters
pa.lambdaH = 10; % penalty coefficient
pa.kap1per2 = 0; % kap1=kap2=1/2 or not? "1"=yes, "0"=no

%% Ghost penalty
%------------------------------------
pa.useGP = 0; % wanna use ghost penalty term?
% There is some bugs (fix later)
pa.gam1 = 1; % parameter for 1st term
pa.gam2 = 1; % parameter for 2nd term

if findCR==1 % wanna find the convergence rate
    nStep = 4; % number of interations
    pa.reguMesh =0;
    wannaPlot = 0; % DON'T CHANGE!
    useDOF = 0; % use dofs to find convergence rate? (default = h)
else % should be used only for plotting
    nStep = 1; % shouldn't change
    pa.reguMesh = 1; % use regular mesh or not?
    nSeg = 15; % ONLY FOR nStep=1;
    wannaPlot = 1; % wanna plot or not? (JUST FOR nStep=1)
    % Go to the end of this file to change what you wanna plot
end

%% settings for each model
%------------------------------------
switch mdl
    case 1 % sinha
        model=sinha;
        pa.kk1 = 1; pa.kk2 = 0.5; % diffusion coefficients (kk1 must be = 2kk2)
    case 2 % barrau (p.44)
        model=barrau;
        pa.kk1 = 1; pa.kk2 = 1000; % for barrau case
        pa.xi = 0.75; % only used for barrau's model, between (-1,1)
    case 3 % barrau (p.47)
        model=barrau_r;
        pa.kk1 = 1; pa.kk2 = 1000; % for barrau case
        pa.xi = 0.75; % only used for barrau's model, between (0,1)
end
GeoDom = model.domain(); % domain


%% ========================================================
% VARIABLES
% =========================================================
numSeg = zeros(nStep,1); % number of segments
hEdgeMax = zeros(nStep,1); % maximum edge size (edge of domain)
phi = cell(nStep,1); % level set function for plotting
% err.L2, err.H1, err.L2G: error in L2, H1 and seminorm.
% err.ENorm: error in ENorm (norm given in Barrau's thesis, p.33)
err = struct('L2',cell(nStep,1),'H1',cell(nStep,1),'L2G',cell(nStep,1),...
        'ENorm',cell(nStep,1),'NJU',cell(nStep,1),...
        'etaK',cell(nStep,1),'etaS',cell(nStep,1),'zetaS',cell(nStep,1));
% sol.ex : exact solution
% sol.Vh: numerical solution in standard FEM
% sol.num: numerical solution in NXFEM
sol = struct('ex',cell(nStep,1),'Vh',cell(nStep,1),'num',cell(nStep,1),'errSolVh',cell(nStep,1));
% msh.hTmax: maximum of all diameters
% msh.ndof: number of dofs.
% msh.p, msh.t, msh.e: points, triangles, edges
msh = struct('p',cell(nStep,1),'t',cell(nStep,1),'e',cell(nStep,1),...
      'hTmax',cell(nStep,1),'ndof',cell(nStep,1),'newNodes',cell(nStep,1));
% just for testing
% test = struct('ii',cell(nStep,1),'jj',cell(nStep,1),'vv',cell(nStep,1));


%% ========================================================
% MAIN PROCESS
% =========================================================
for z=1:nStep
    if nStep>1
%         numSeg(z) = 2^z*10+1; % 21,41,81,161
        numSeg(z) = 2^(z-1)*10+1; % 11,21,41,81
%         numSeg = [21,41,81,121];
    else % nStep=1
        numSeg(z) = nSeg;
    end
    hEdgeMax(z) = 2/numSeg(z); % maximum edge size (edge of domain)
    [msh(z),err(z),sol(z),phi{z}] = main_eachStep(hEdgeMax(z),GeoDom,pa,model);
end


tx=zeros(1,nStep); ty=zeros(4,nStep);
%% ========================================================
% CONVERGENCE RATE
% display convergence rate
% =========================================================
if nStep>1
    for i=1:nStep
       tx(i)=msh(i).hTmax;
       ty(1,i)=err(i).L2;
       ty(2,i)=err(i).H1;
       ty(3,i)=err(i).L2G;
       ty(4,i)=err(i).ENorm;
    end
    tmp = polyfit(log(tx),log(ty(1,:)),1);
    order.L2 = tmp(1);
    tmp = polyfit(log(tx),log(ty(2,:)),1);
    order.H1 = tmp(1);
    tmp = polyfit(log(tx),log(ty(3,:)),1);
    order.L2G = tmp(1);
    tmp = polyfit(log(tx),log(ty(4,:)),1);
    order.ENorm = tmp(1);
    
    % plot results
    plot(log(tx),log(ty(1,:)),'-r.',...
        log(tx),log(ty(4,:)),'-b.');
    legend('L2','ENorm');
    xlabel('log(h)'); 
    ylabel('log(error)');
    
    % display results
    fprintf('regular: %d\n',pa.reguMesh);
    fprintf('kap1per2: %d\n',pa.kap1per2);
    fprintf('numSeg: %0.0f\n',numSeg(nStep));
    fprintf('DOFs: %0.0f\n',msh(nStep).ndof);
    fprintf('L2 order: %0.5f\n',order.L2);
    fprintf('L2Grad order: %0.5f\n',order.L2G);
    fprintf('H1 order: %0.5f\n',order.H1);
    fprintf('ENorm order: %0.5f\n',order.ENorm);
    
%     z = nStep-1;
%     if ~useDOF
%         % via msh.hTmax
%         order.L2 = log(err(z+1).L2/err(z).L2)/log(msh(z+1).hTmax/msh(z).hTmax);
%         order.H1 = log(err(z+1).H1/err(z).H1)/log(msh(z+1).hTmax/msh(z).hTmax);
%         order.L2G = log(err(z+1).L2G/err(z).L2G)/log(msh(z+1).hTmax/msh(z).hTmax);
%         order.ENorm = log(err(z+1).ENorm/err(z).ENorm)/log(msh(z+1).hTmax/msh(z).hTmax);
%     else
%         % via msh.ndof
%         order.L2 = 2*log(err(z+1).L2/err(z).L2)/log(msh(z).ndof/msh(z+1).ndof);
%         order.H1 = 2*log(err(z+1).H1/err(z).H1)/log(msh(z).ndof/msh(z+1).ndof);
%         order.L2G = 2*log(err(z+1).L2G/err(z).L2G)/log(msh(z).ndof/msh(z+1).ndof);
%         order.ENorm = 2*log(err(z+1).ENorm/err(z).ENorm)/log(msh(z).ndof/msh(z+1).ndof);
%     end
% 
%     % display results
%     fprintf('regular: %d\n',pa.reguMesh);
%     fprintf('kap1per2: %d\n',pa.kap1per2);
%     fprintf('numSeg: %0.0f\n',numSeg(z+1));
%     fprintf('DOFs: %0.0f\n',msh(z+1).ndof);
%     fprintf('L2 order: %0.5f\n',order.L2);
%     fprintf('L2Grad order: %0.5f\n',order.L2G);
%     fprintf('H1 order: %0.5f\n',order.H1);
%     fprintf('ENorm order: %0.5f\n',order.ENorm);
else
    fprintf('numSeg: %0.0f\n',numSeg(nStep));
    fprintf('DOFs: %0.0f\n',msh(nStep).ndof);
    fprintf('L2 err: %0.8f\n',err(nStep).L2);
    fprintf('L2G err: %0.8f\n',err(nStep).L2G);
    fprintf('H1 err: %0.8f\n',err(nStep).H1);
    fprintf('ENorm err: %0.8f\n',err(nStep).ENorm);
    fprintf('etaK: %0.5f\n',err(nStep).etaK);
    fprintf('etaS: %0.5f\n',err(nStep).etaS);
    fprintf('zetaS: %0.5f\n',err(nStep).zetaS);
    fprintf('NJU: %0.5f\n',err(nStep).NJU);
end



%% ========================================================
% PLOTTING
% =========================================================
if wannaPlot==1
    D2=2; D3=3; D2D3=1; % which dimension? (D1D2 for both cases)
    wM=1; wtM=0; % plot with mesh or without mesh?
    On = 'on'; Off = 'off';
    % MESH --------------------------------------
    mshs{1} = 1; % wanna plot mesh? (1 or 0)
    mshs{2} = On; % node labels? (On or Off)
    mshs{3} = On; % element label? (On or Off)
    % LEVEL SET FUNCTION ------------------------
    lsf{2} = 0; % wanna plot level set function? (1 or 0)
    lsf{1} = phi{1}; % level set function
    lsf{3} = D2; % dimension? (D2, D3 or D1)
    lsf{4} = wtM; % plot with mesh? (wM or wtM)
    % EXACT SOLUTION ----------------------------
    eS{2} = 1; % wanna plot exact solution? (1 or 0)
    eS{1} = sol(1).ex; % exact solution
    eS{3} = D3; % dimension? (D2, D3 or D1)
    eS{4} = wM; % plot with mesh? (wM or wtM)
    % NUMERICAL SOLUTION ------------------------
    nS{2} = 1; % wanna plot numerical solution? (1 or 0)
    nS{1} = sol(1).Vh; % numerical solution (in standard FEM)
%     nS{1} = sol(1).errSolVh; % error |u-uex| (in standard FEM)
    nS{3} = D3; % dimension? (D2, D3 or D1)
    nS{4} = wM; % plot with mesh? (wM or wtM)
    % DO THE PLOT -------------------------------------
    plotNXFEM(nS,eS,mshs,lsf,msh);
end
