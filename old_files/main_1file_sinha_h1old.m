%% ========================================================
% INTRODUCTION
% =========================================================
% NXFEM based on Hansbo's paper.
% This file is used for test/plotting the solution/compute condition number
% Test case given in Sinha's paper
% Omg=[0,2]x[0,1], Gam={1}x[0,1]
%   -nabla\cdot(k\nabla u) = f in Omg_i
%   [u]=[k\nabla_n u]=0 on Gam
%   u=0 on \part\Omg
% Exact solution :  given in defExSol.m
% k1 = 2k2
% =========================================================


%% ========================================================
% GLOBAL PARAMETERS
% =========================================================
global kk1 kk2 lambdaH degP1D degP2D degN kap1per2 % input parameters



%% ========================================================
% INPUT
% =========================================================
degP1D = 3; % Gaussian quadrature points in 1D (for polinomial functions)
degP2D = 4; % Gaussian quadrature points in 2D (for polinomial functions)
degN = 8; % Gaussian quadrature points in 2D (for non-polynomial functions)
% degree-#OfPoints : 1-1, 2-3, 3-4, 4-6, 5-7, 6-12, 
%                    7-13, 8-16, 9-19, 10-25, 11-27, 12-33

xDomVal = [0 2 2 0]; % x values of points constructing Omega
yDomVal = [0 0 1 1]; % corresponding y value
RectDom = [3,4,xDomVal,yDomVal]'; % rectangular domain "3" with "4" sides
GeoDom = decsg(RectDom);

numSeg = 40; % number of segments
hEdgeMax = 2/numSeg; % maximum edge size (edge of domain)

kk1 = 1; kk2 = 0.5; % diffusion coefficients
lambdaH = 1e2 ; 
kap1per2 = 0; % true then kap1=kap2=1/2

% ghost penalty
global gam1 gam2
gam1 = 1; 
gam2 = 1;



%% ========================================================
% PLOTTING
% =========================================================
D2=2; D3=3; D2D3=1; % which dimension? (D1D2 for both cases)
wM=1; wtM=0; % plot with mesh or without mesh?
On = 'on'; Off = 'off';
% MESH --------------------------------------
mSh{1} = 0; % wanna plot mesh? (1 or 0)
mSh{2} = Off; % node labels? (On or Off)
mSh{3} = Off; % element label? (On or Off)
% LEVEL SET FUNCTION ------------------------
lsf{2} = 0; % wanna plot level set function? (1 or 0)
lsf{3} = D3; % dimension? (D2, D3 or D1)
lsf{4} = wtM; % plot with mesh? (wM or wtM)
% EXACT SOLUTION ----------------------------
eS{2} = 1; % wanna plot exact solution? (1 or 0)
eS{3} = D3; % dimension? (D2, D3 or D1)
eS{4} = wtM; % plot with mesh? (wM or wtM)
% NUMERICAL SOLUTION ------------------------
nS{2} = 1; % wanna plot numerical solution? (1 or 0)
nS{3} = D3; % dimension? (D2, D3 or D1)
nS{4} = wtM; % plot with mesh? (wM or wtM)



%% ========================================================
% MESH
% =========================================================
% get information of the mesh
global points edges triangles
% [points,edges,triangles] = getMeshFromFF('ffmeshSinha.msh'); % use ff++
% [points,edges,triangles] = initmesh(GeoDom,'hmax',hEdgeMax); % use matlab
[points,edges,triangles] = poimesh(GeoDom,2*numSeg+1,numSeg); % regular mesh
x = points(1,:); % x-coordinate of points
y = points(2,:); % y-coordinate of points
hT = getDiam(points,triangles); % diameter (longest side) of each triangle
hTmax = max(hT); % maximum of all diameters

% level set function (interface)
phi = defPhi(x,y); % row-array



%% ========================================================
% GET INFORMATION
% =========================================================
% TRIANGLES
[CTs,omg1NCTs,omg2NCTs,areaNCTs1,areaNCTs2,areaCTs] = getTriangles(phi);

% ON CUT TRIANGLES
[idxEachCTs,idxCTs,typeCTs,areaChildCTs,iPs,uNormalsCT]...
            = getInfoCTs(CTs,areaCTs,phi);

% NODES
global newNodes
nNewNodes = idxCTs{4}; % number of new nodes (nodes around the interface)
nNodes = size(points,2); % number of standard nodes

% vector contaning new numbering of nodes around interface
newNodes = getNewNodes(idxCTs{5},nNodes); % column-array
% nodes in CTs, row-array
nodesInCTs = idxCTs{5};
% nodes of NCTs (nodes in Omg1 and nodes on Gam), column-array
nodesOmg1Gam = unique(omg1NCTs(1:3,:)); % column-array
% nodes of NCTs (nodes in Omg2 and nodes on Gam), column-array
nodesOmg2Gam = unique(omg2NCTs(1:3,:)); % column-array
% nodes in Omg2 but not in CTs (colum-array)
% there may be still some nodes on the interface
nodesOmg2NotInCTs = setdiff(nodesOmg2Gam,nodesInCTs);
% nodes in Omg2 (not on interface) of cut triangles
nodesCTsInOmg2NotGam = idxCTs{3}; % row-array
 % nodes in Omg1 (not on interface CTs)
nodesOmg1NotGamCTs = setdiff(nodesOmg1Gam,idxCTs{1});
nodesOnGamOfCTs = idxCTs{1}; % nodes on interface of CTs   
% nodes of CTs in Omg2 and on Gamma
nodesCTsInOmg2OnGam = [idxCTs{3},idxCTs{1}]; % row-array

% just to be used for "another way to find H1"
%----------------------------------------------
% nodes in Omg1 and nodes in CTs (colum-array)
nodesOmg1AndInCTs = unique([nodesOmg1Gam;nodesInCTs']);
% nodes in Omg2 and nodes in CTs (colum-array)
nodesOmg2AndInCTs = unique([nodesOmg2Gam;nodesInCTs']);



%% ========================================================
% CONTROL PARAMETERS
% depend on mesh
% =========================================================
global lambda kap1 kap2
[kap1,kap2] = getKappa(areaChildCTs); % follows Becker's guide
lambda = getLambda(areaChildCTs,CTs,iPs,hT,hTmax); % follows Burman's guide



%% ========================================================
% STIFFNESS MATRIX
% on all nodes including nodes on the boundary
% =========================================================
A = getGlobalMatrix(omg1NCTs,omg2NCTs,CTs,phi,hTmax, ...
                 areaChildCTs,areaCTs,iPs,uNormalsCT,nodesCTsInOmg2OnGam);



%% ========================================================
% LOAD VECTOR
% the right hand side
% =========================================================
F = getLoad(omg1NCTs,omg2NCTs,areaNCTs1,areaNCTs2,...
                    CTs,iPs,areaCTs,typeCTs,idxEachCTs,idxCTs);



%% ========================================================
% BOUNDARY CONDITIONS
% =========================================================
bNodes = unique([edges(1,:) edges(2,:)]); % boundary nodes,
% boundary nodes around interface
bNodesInCTs = intersect(bNodes,nodesInCTs); 
% NEW boundary nodes around interface
bNodesInCTs = newNodes(bNodesInCTs); 
bNodes = [bNodes,bNodesInCTs']; % update boundary nodes

iNodes = setdiff(1:nNodes,bNodes);
% inner nodes around interface
iNodesInCTs = intersect(iNodes,nodesInCTs); 
% NEW boundary nodes around interface
iNodesInCTs = newNodes(iNodesInCTs);
iNodes = [iNodes,iNodesInCTs']; % update inner nodes



%% ========================================================
% SEEK NUMERICAL SOLUTION
% in NXFEM space
% =========================================================
sol = zeros(nNodes+nNewNodes,1); % column-array
sol(bNodes) = 0; % homogeneous Dirichlet BC
% LU factorization
sol(iNodes) = A(iNodes,iNodes)\F(iNodes); % don't care nodes on boundary

% GMRES factorization
% sol(iNodes) = gmres(A(iNodes,iNodes),F(iNodes));

BB = A(iNodes,iNodes);
condB = condest(BB);
fprintf('condB: %0.5f\n',condB);


%% ========================================================
% NUMERICAL SOLUTION 
% in STANDARD FEM (just for plotting)
% sol4plot_i = sol_i for i in nodesOmg1NotOnGam or nodesOmg2NotAroundGam
% sol4plot_i = sol_k(i) for i in nodesOmg2CT
% sol4plot_i = sol_i+sol_k(i) for i in nodesOnGam
% =========================================================
sol4plot(nodesOmg1NotGamCTs) = sol(nodesOmg1NotGamCTs);
sol4plot(nodesOmg2NotInCTs) = sol(nodesOmg2NotInCTs);
sol4plot(nodesCTsInOmg2NotGam) = sol(newNodes(nodesCTsInOmg2NotGam));
sol4plot(nodesOnGamOfCTs) = sol(newNodes(nodesOnGamOfCTs))...
                                + sol(nodesOnGamOfCTs);
sol4plot = sol4plot'; % transform to column-array for pdesurf



%% ========================================================
% EXACT SOLUTION
% in standard FEM space (just for plotting)
% exSol_i = exSol(x_i)
% =========================================================
exSol = zeros(nNodes,1); % column-array
exSol(nodesOmg1Gam) = defExSol(points(1,nodesOmg1Gam),...
        points(2,nodesOmg1Gam),1);
exSol(nodesOmg2NotInCTs) = defExSol(points(1,nodesOmg2NotInCTs),...
        points(2,nodesOmg2NotInCTs),2);
exSol(nodesCTsInOmg2NotGam) = defExSol(points(1,nodesCTsInOmg2NotGam),...
        points(2,nodesCTsInOmg2NotGam),2);



%% ========================================================
% EXACT SOLUTION
% in NXFEM (for finding |exSol-sol|)
% exSolNX_i = exSol_i for i is node of mesh
% exSolNX_k(i) = exSol_i for i in nodesAroundGam
% =========================================================
exSolNX = zeros(nNewNodes,1); % column-array
allNodes = 1:nNodes;
exSolNX(allNodes) = exSol(allNodes);
exSolNX(newNodes(nodesInCTs)) = exSol(nodesInCTs);



%% ========================================================
% A PRIORI ERROR ESTIMATES
% error between numerical & exact solutions (column-array)
% =========================================================
err = exSolNX - sol; % column-array

% ---------------------------------------------------------
% L2
% ---------------------------------------------------------
mL2 = getMatrixL2(omg1NCTs,omg2NCTs,areaNCTs1,areaNCTs2,...
                iPs,CTs,typeCTs,idxEachCTs,areaCTs,nodesCTsInOmg2OnGam);
errL2square = err' * mL2 * err;
errL2 = sqrt(errL2square);

% ---------------------------------------------------------
% H1 : We need to find ||Grad e||_L2 first
% ---------------------------------------------------------
mL2gradOld = getMatrixL2Grad(omg1NCTs,omg2NCTs,nodesCTsInOmg2OnGam,...
                    CTs,areaChildCTs,areaCTs);
errL2grad = err' * mL2gradOld * err;
errH1Old = errL2square + errL2grad;
errH1Old = sqrt(errH1Old);

% ---------------------------------------------------------
% H1
% another way to find H1
% ---------------------------------------------------------
% Get numerical solution in Omg1
nMaxOmg1 = max(nodesOmg1AndInCTs);
err1 = zeros(nMaxOmg1,1); % column-array
err1(nodesOmg1AndInCTs) = err(nodesOmg1AndInCTs);
% Get numerical solution on Omg2
nMaxOmg2 = max(nodesOmg2AndInCTs);
err2 = zeros(nMaxOmg2,1); % column-array
err2(nodesOmg2NotInCTs) = err(nodesOmg2NotInCTs);
err2(nodesInCTs) = err(newNodes(nodesInCTs));

[mH11,mH12] = getMatrixH1(omg1NCTs,omg2NCTs,areaNCTs1,...
               areaNCTs2,CTs,iPs,typeCTs,areaChildCTs,areaCTs,idxEachCTs);
errH1 = err1'*mH11*err1 + err2'*mH12*err2;
errH1 = sqrt(errH1);

% print out
fprintf('L2 err: %0.5f\n',errL2);
fprintf('H1 err: %0.5f\n',errH1);
fprintf('H1old err: %0.5f\n',errH1old);


%% ========================================================
% PLOTTING
% =========================================================
% LEVEL SET FUNCTION ------------------------
lsf{1} = phi; % level set function
% EXACT SOLUTION ----------------------------
eS{1} = exSol; % exact solution
% NUMERICAL SOLUTION ------------------------
nS{1} = sol4plot; % numerical solution (in standard FEM)
% DO THE PLOT -------------------------------------
plotNXFEM(nS,eS,mSh,lsf);
