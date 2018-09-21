%% ========================================================
% INTRODUCTION
% =========================================================
% NXFEM based on Hansbo's paper.
% This file is used to find the convergece rate
% Test case given in Sinha's paper
% Omg=[0,2]x[0,1], Gam={1}x[0,1]
%   -nabla\cdot(k\nabla u) = f in Omg_i
%   [u]=[k\nabla_n u]=0 on Gam
%   u=0 on \part\Omg
% Exact solution :  given in defExSol.m
% k1 = 2k2 (obligatory)
% =========================================================


%% ========================================================
% FIXED PARAMERTERS
% shouldn't change
% =========================================================
global degP1D degP2D degN
degP1D = 3; % Gaussian quadrature points in 1D (for polinomial functions)
degP2D = 4; % Gaussian quadrature points in 2D (for polinomial functions)
degN = 8; % Gaussian quadrature points in 2D (for non-polynomial functions)
% degree-#OfPoints : 1-1, 2-3, 3-4, 4-6, 5-7, 6-12, 
%                    7-13, 8-16, 9-19, 10-25, 11-27, 12-33



%% ========================================================
% INPUT REGION
% =========================================================
global kk1 kk2 lambdaH kap1per2 xi% input parameters
global model % indicate the model used for the whole project
% available models: sinha, barrau, carraro 
model=sinha;

GeoDom = model.domain(); % domain

findCR = 0; % wanna find convergence rate? 1 or 0
lambdaH = 1; % penalty coefficient
kap1per2 = 0; % kap1=kap2=1/2 or not? "1"=yes, "0"=no
kk1 = 1; kk2 = 0.5; % diffusion coefficients (kk1 must be = 2kk2)
xi = 1; % not used for sinha's case

global gam1 gam2 useGP
useGP = 0; % wanna use ghost penalty term?
gam1 = 1;
gam2 = 1;

if findCR==1 %
   nStep = 2; % DON'T CHANGE! 
   reguMesh = 1; % DON'T CHANGE!
   wannaPlot = 0; % DON'T CHANGE!
else % should be used only for plotting
    nStep = 1; % 1 or 2
    nSeg = 41; % ONLY FOR nStep=1;
    wannaPlot = 1; % wanna plot or not? (JUST FOR nStep=1)
    reguMesh = 1; % use regular mesh or not? (JUST FOR nStep=1)
end


%% ========================================================
% VARIABLES
% =========================================================
numSeg = zeros(nStep,1); % number of segments
hEdgeMax = zeros(nStep,1); % maximum edge size (edge of domain)
hTmax = zeros(nStep,1); % max diam of all triangles
errL2  = zeros(nStep,1); % err in L2, column array
errH1  = zeros(nStep,1); % err in H1, column array
errL2G  = zeros(nStep,1); % grad err in L2, column array
nDOFs  = zeros(nStep,1); % number of DOFs, column array
sol4plot = cell(nStep,1); % numerical solution for plotting
exSol = cell(nStep,1); % exact solution for plotting
phi = cell(nStep,1); % level set function for plotting



%% ========================================================
% MAIN PROCESS
% =========================================================
for z=1:nStep
    if nStep>1
        numSeg(z) = 40*z+1; % just for test
%         numSeg(z) = 2^(z+2)*10+1; % number of segments (57,113)
    else % nStep=1
        numSeg(z) = nSeg;
    end
    hEdgeMax(z) = 2/numSeg(z); % maximum edge size (edge of domain)
    [nDOFs(z),hTmax(z),errL2(z),errH1(z),sol4plot{z},exSol{z},phi{z}]...
            = getErrEachStep(hEdgeMax(z),GeoDom,reguMesh);
end



%% ========================================================
% CONVERGENCE RATE
% display convergence rate
% =========================================================
if nStep>1
    z = nStep-1;
    % via hTmax
    orderL2 = log(errL2(z+1)/errL2(z))/log(hTmax(z+1)/hTmax(z)); 
    orderH1 = log(errH1(z+1)/errH1(z))/log(hTmax(z+1)/hTmax(z)); 
    orderL2G = log(errL2G(z+1)/errL2G(z))/log(hTmax(z+1)/hTmax(z)); 
    
    % via nDOFs
%     orderL2 = 2*log(errL2(z+1)/errL2(z))/log(nDOFs(z)/nDOFs(z+1)); 
%     orderH1 = 2*log(errH1(z+1)/errH1(z))/log(nDOFs(z)/nDOFs(z+1)); 
%     orderL2G = 2*log(errL2G(z+1)/errL2G(z))/log(nDOFs(z)/nDOFs(z+1)); 
    
    % display results
    fprintf('lambdaH: %0.0f\n',lambdaH);
    fprintf('DOFs: %0.0f\n',nDOFs(z+1));
    fprintf('L2 order: %0.5f\n',orderL2);
%     fprintf('L2Grad order: %0.5f\n',orderL2G);
    fprintf('H1 order: %0.5f\n',orderH1);
else
    fprintf('L2 err: %0.5f\n',errL2(nStep));
    fprintf('H1 err: %0.5f\n',errH1(nStep));
end



%% ========================================================
% PLOTTING
% =========================================================
if wannaPlot==1
    D2=2; D3=3; D2D3=1; % which dimension? (D1D2 for both cases)
    wM=1; wtM=0; % plot with mesh or without mesh?
    On = 'on'; Off = 'off';
    % MESH --------------------------------------
    mSh{1} = 0; % wanna plot mesh? (1 or 0)
    mSh{2} = Off; % node labels? (On or Off)
    mSh{3} = Off; % element label? (On or Off)
    % LEVEL SET FUNCTION ------------------------
    lsf{2} = 0; % wanna plot level set function? (1 or 0)
    lsf{1} = phi{1}; % level set function
    lsf{3} = D2; % dimension? (D2, D3 or D1)
    lsf{4} = wtM; % plot with mesh? (wM or wtM)
    % EXACT SOLUTION ----------------------------
    eS{2} = 1; % wanna plot exact solution? (1 or 0)
    eS{1} = exSol{1}; % exact solution
    eS{3} = D3; % dimension? (D2, D3 or D1)
    eS{4} = wtM; % plot with mesh? (wM or wtM)
    % NUMERICAL SOLUTION ------------------------
    nS{2} = 1; % wanna plot numerical solution? (1 or 0)
    nS{1} = sol4plot{1}; % numerical solution (in standard FEM)
    nS{3} = D3; % dimension? (D2, D3 or D1)
    nS{4} = wtM; % plot with mesh? (wM or wtM)
    % DO THE PLOT -------------------------------------
    plotNXFEM(nS,eS,mSh,lsf); 
end