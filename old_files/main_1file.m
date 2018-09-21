%% ========================================================
% INTRODUCTION
% =========================================================
% NXFEM based on Hansbo's paper.
% This file is used for testing, ploting,...
% (not used for finding convergence rate, pls use main.m instead)
% This file is used for both Sinha's models and Barrau's model
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

% available models: sinha, barrau 
% model=sinha;
model=barrau;

GeoDom = model.domain(); % domain
numSeg = 41; % number of segments
hEdgeMax = 2/numSeg; % maximum edge size (edge of domain)

lambdaH = 10; % penalty coefficient
kap1per2 = 0; % kap1=kap2=1/2 or not? "1"=yes, "0"=no
% kk1 = 1; kk2 = 0.5; % diffusion coefficients (for sinha's case)
kk1 = 1; kk2 = 1000; % for barau's case;
xi = 0.71; % only used for barrau's model
reguMesh = 0; % use regular mesh or not? (JUST FOR nStep=1)

% Ghost penalty
global gam1 gam2 useGP
useGP = 0; % wanna use ghost penalty term?
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
mSh{2} = On; % node labels? (On or Off)
mSh{3} = On; % element label? (On or Off)
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
if reguMesh
    [points,edges,triangles] = poimesh(GeoDom,2*numSeg+1,numSeg); % regular mesh
else
    [points,edges,triangles] = initmesh(GeoDom,'hmax',hEdgeMax); % use matlab
end
x = points(1,:); % x-coordinate of points
y = points(2,:); % y-coordinate of points
hT = getDiam(points,triangles); % diameter (longest side) of each triangle
hTmax = max(hT); % maximum of all diameters

% level set function (interface)
defPhi = model.defPhi;
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
% EXACT SOLUTION
% in standard FEM space (just for plotting)
% exSol_i = exSol(x_i)
% =========================================================
defExSol = model.defExSol;
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
% SEEK NUMERICAL SOLUTION
% in NXFEM space
% =========================================================
sol = zeros(nNodes+nNewNodes,1); % column-array
typeBC = model.bc(); % get type of BCs
switch typeBC
    case 1 % u=o on whole boundary
        sol(bNodes) = 0;
    case 2 % u=uex on whole boundary
        sol(bNodes) = exSolNX(bNodes);
end
% modification of F
F = F - A*sol;
% LU factorization
sol(iNodes) = A(iNodes,iNodes)\F(iNodes); % don't care nodes on boundary
% GMRES factorization
% sol(iNodes) = gmres(A(iNodes,iNodes),F(iNodes));



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
mL2grad = getMatrixL2Grad(omg1NCTs,omg2NCTs,nodesCTsInOmg2OnGam,...
                    CTs,areaChildCTs,areaCTs);
errL2grad = err' * mL2grad * err;
errH1 = errL2square + errL2grad;
errH1 = sqrt(errH1);

% print out
fprintf('L2 err: %0.5f\n',errL2);
fprintf('H1 err: %0.5f\n',errH1);



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



%% ========================================================
% POSTERIORI ESTIMATES
% =========================================================

% f in standard FEM 
% ------------------------
% (in this model, f is the same for both subdomains but we still use these
% lines of code)
defF = model.defF;
fSD = zeros(nNodes,1); % column-array
fSD(nodesOmg1Gam) = defF(points(1,nodesOmg1Gam),...
        points(2,nodesOmg1Gam),1);
fSD(nodesOmg2NotInCTs) = defF(points(1,nodesOmg2NotInCTs),...
        points(2,nodesOmg2NotInCTs),2);
fSD(nodesCTsInOmg2NotGam) = defF(points(1,nodesCTsInOmg2NotGam),...
        points(2,nodesCTsInOmg2NotGam),2);

% f in NXFEM 
% ------------------------
fNX = zeros(nNewNodes,1); % column-array
fNX(allNodes) = fSD(allNodes);
fNX(newNodes(nodesInCTs)) = fSD(nodesInCTs);

% etaK
% -----------------------------
metaK = getMatEstEtaK(omg1NCTs,omg2NCTs,areaNCTs1,areaNCTs2,...
             iPs,CTs,typeCTs,idxEachCTs,areaCTs,nodesCTsInOmg2OnGam,hT);
etaKSquare = fNX'*metaK*fNX;
etaKSquare