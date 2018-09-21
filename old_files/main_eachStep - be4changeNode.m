function [msh,err,sol,pplot,phi] = main_eachStep(hEdgeMax,GeoDom,pa,model)
% Find L2 err and H1 err for each step
% This funciton is used for finding convergence rate in main_CR.m
% Input: hEdgeMax, GeoDown to create a mesh
% Output: - error in L2 norm and H1 norm wrt this step
%         - some testing outputs
%         - some useful info about the mesh and solution


%% ========================================================
% MESH
% =========================================================
% get information of the mesh
if ~pa.reguMesh
    [points,edges,triangles] = initmesh(GeoDom,'hmax',hEdgeMax); % irregular
else
    numSeg = 2/hEdgeMax;
%     [points,edges,triangles] = poimesh(GeoDom,2*numSeg+1,numSeg); % regular
    [points,edges,triangles] = poimesh(GeoDom,numSeg,numSeg); % regular
end
% save to msh
msh.p = points; msh.t = triangles; msh.e = edges;
x = points(1,:); % x-coordinate of points
y = points(2,:); % y-coordinate of points
% diameter (longest side) of each triangle: 1 x nTs
hT = getDiam(msh); 
msh.hTmax = max(hT); % maximum of all diameters
tol=pa.tol;
% level set function (interface)
phi = model.defPhi(x,y,pa); % 1 x number of points (row array)
phi(abs(phi)<tol)=0; % find values of phi which are very small (~0) and set them to 0
% exact solution function
defExSol = model.defExSol;



%% ========================================================
% GET INFORMATION OF TRIANGLES
% =========================================================

%-------------------------------------------------------------------
% TRIANGLES
%--------------------------------------------
[CTs,NCTs1,NCTs2] = getTriangles(phi,msh,pa);

%-------------------------------------------------------------------
% ON CUT TRIANGLES
%-------------------------------------------------------------------
[nodeCTs,typeCTs,areaChildCTs,iPs,uNVCTs,areaCTs]...
            = getInfoCTs(CTs,phi,msh,pa);
        
%-------------------------------------------------------------------
% Find small-cut triangles (idx in the OLD CTs)
%-------------------------------------------------------------------
if pa.smallCut
    [tri2del,t2Omg1,t2Omg2] = findSmallPhi(msh,CTs,iPs,hT,...
                                    nodeCTs,pa,typeCTs,areaCTs,phi);
    % If there are small-cut triangles, remove them from CTs!!
    if ~isempty(tri2del)
        nCTs = size(CTs,2); % number of OLD cut triangles
        % get NEW not-cut triangles
        if ~isempty(t2Omg1)
            NCTs1 = [NCTs1,CTs(:,t2Omg1)]; % add more triangles to nCTs1
        end
        if ~isempty(t2Omg2)
            NCTs2 = [NCTs2,CTs(:,t2Omg2)]; % add more triangles to nCTs2
        end
        % get NEW cut triangles
        CTs = CTs(:,setdiff(1:nCTs,tri2del));
        % find again all information
        clear nodeCTs typeCTs areaChildCTs iPs uNVCTs areaCTs; % just in case
        [nodeCTs,typeCTs,areaChildCTs,iPs,uNVCTs,areaCTs]...
                = getInfoCTs(CTs,phi,msh,pa);
    end
end

%--------------------------------------------
% area of not cut triangles
%--------------------------------------------
areaNCTs1 = getAreaTris(NCTs1,msh);
areaNCTs2 = getAreaTris(NCTs2,msh);



%% ========================================================
% NODES
% =========================================================
% get info
nodeCTsOnG = nodeCTs.onG;
nodeCTsOmg2 = nodeCTs.Omg2;
nodeCTsall = nodeCTs.all;
nNewNodes = nodeCTs.n; % number of new nodes (nodes around the interface)
nNodes = size(points,2); % number of standard nodes
msh.ndof = nNewNodes + nNodes; % number of dofs
stdNodes = 1:nNodes; % vecotor of standard nodes
nodesOmg1 = stdNodes((phi<0)&(abs(phi)>tol)); % all std nodes inside Omg1
nodesOnGam = stdNodes(abs(phi)<tol); % all std nodes on Gam
nodesOmg2 = stdNodes((phi>0)&(abs(phi)>tol)); % all std nodes inside Omg2
% vector contaning new numbering of nodes around interface
newNodes = getNewNodes(nodeCTsall,nNodes); % column-array
msh.newNodes=newNodes;
% nodes in CTs, row-array
nodesInCTs = nodeCTsall;
% nodes of NCTs (nodes in Omg1 and nodes on Gam), column-array
nodesOmg1Gam = unique(NCTs1(1:3,:)); % column-array, there still nodes on Gam in CTs which not belongs to NCTs
nodesOmg1Gam = union(nodesOmg1Gam,nodeCTsOnG); % column-array
% nodes of NCTs (nodes in Omg2 and nodes on Gam), column-array
nodesOmg2Gam = unique(NCTs2(1:3,:)); % column-array, there still nodes on Gam in CTs which not belongs to NCTs
nodesOmg2Gam = union(nodesOmg2Gam,nodeCTsOnG); % column-array
% nodes in Omg2 but not in CTs (colum-array)
% there may be still some nodes on the interface
nodesOmg2NotInCTs = setdiff(nodesOmg2Gam,nodesInCTs);
% nodes in Omg2 (not on interface) of cut triangles
nodesCTsInOmg2NotGam = nodeCTsOmg2; % row-array
 % nodes in Omg1 (not on interface CTs)
nodesOmg1NotGamCTs = setdiff(nodesOmg1Gam,nodeCTsOnG);
nodesOnGamOfCTs = nodeCTsOnG; % nodes on interface of CTs   
% nodes in CTs on Omg2 and on Gamma
nodesCTsInOmg2OnGam = [nodeCTsOmg2,nodeCTsOnG]; % row-array



%% ========================================================
% CONTROL PARAMETERS
% depend on mesh
% =========================================================
% cp = changing parameters
[cp.kap1,cp.kap2] = getKappa(areaChildCTs,pa); % follows Becker's guide
% cp.lambda = getLambda(areaChildCTs,CTs,iPs,hT,msh.hTmax,pa);
cp.lambda = getLambda(areaChildCTs,CTs,iPs,hT,msh.hTmax,pa,areaCTs);
cp.kk1 = pa.kk1; cp.kk2 = pa.kk2;


%% ========================================================
% STIFFNESS MATRIX
% on all nodes including nodes on the boundary
% =========================================================
A = getGlobalMatrix(NCTs1,NCTs2,CTs,phi,hT,...
     areaChildCTs,areaCTs,iPs,uNVCTs,nodesCTsInOmg2OnGam,msh,pa,cp);

% A = getGMtwLinda(NCTs1,NCTs2,CTs,phi,hT,...
%      areaChildCTs,areaCTs,iPs,uNVCTs,nodesCTsInOmg2OnGam,msh,pa,cp);

%% ========================================================
% LOAD VECTOR
% the right hand side
% =========================================================
defF = model.defF;
F = getLoad(NCTs1,NCTs2,areaNCTs1,areaNCTs2,...
                    CTs,iPs,areaCTs,typeCTs,nodeCTs,msh,pa,defF);


 
%% ========================================================
% BOUNDARY CONDITIONS
% =========================================================
bNodesStd = unique([edges(1,:) edges(2,:)]); % standard boundary nodes,
% boundary nodes around interface
bNodesInCTs = intersect(bNodesStd,nodesInCTs);
bNodeInt = bNodesInCTs;
% NEW boundary nodes around interface
bNodesInCTs = newNodes(bNodesInCTs); 
bNodeInt = [bNodeInt,bNodesInCTs']; % std and new nodes on boundary cap interface
bNodes = [bNodesStd,bNodesInCTs']; % all boundary nodes (std + new)
iNodesStd = setdiff(1:nNodes,bNodesStd); % standard inner nodes
% inner nodes around interface
iNodesInCTs = intersect(iNodesStd,nodesInCTs); 
% NEW boundary nodes around interface
iNodesInCTs = newNodes(iNodesInCTs);
iNodes = [iNodesStd,iNodesInCTs']; % all inner nodes (std + new)



%% ========================================================
% EXACT SOLUTION
% in standard FEM space (just for plotting)
% exSol_i = exSol(x_i)
% =========================================================
exSol = zeros(nNodes,1); % column-array
exSol(nodesOmg1) = defExSol(points(1,nodesOmg1),points(2,nodesOmg1),1,pa); % nodes inside Omg1
exSol(nodesOnGam) = defExSol(points(1,nodesOnGam),points(2,nodesOnGam),1,pa); % nodes on Gam
exSol(nodesOmg2) = defExSol(points(1,nodesOmg2),points(2,nodesOmg2),2,pa); % nodes inside Omg2
% add to var
sol.ex = exSol;



%% ========================================================
% EXACT SOLUTION
% in NXFEM (for finding |exSol-sol|)
% exSolNX_i = exSol_i for i is node of mesh
% exSolNX_k(i) = exSol_i for i in nodesAroundGam
% =========================================================
exSolNX = zeros(nNodes+nNewNodes,1); % column-array
exSolNX(stdNodes) = exSol(stdNodes);
exSolNX(newNodes(nodesInCTs)) = exSol(nodesInCTs);



%% ========================================================
% PRECONDITIONER
% nodes outside CTs region
%% ========================================================
useit2=0;
if useit2
    nodesOutsideCTs = setdiff(1:msh.ndof,[nodesInCTs,nNodes:nNodes+nNewNodes]); 
    P = getMatPrecond(nodesOutsideCTs,NCTs1,NCTs2,...
            nodesCTsInOmg2OnGam,CTs,pa,areaChildCTs,areaCTs,msh);
end % useit2



%% ========================================================
% SEEK NUMERICAL SOLUTION
% in NXFEM space
% =========================================================
numSol = zeros(nNodes+nNewNodes,1); % column-array
typeBC = model.bc(); % get type of BCs
switch typeBC
    case 1 % u=o on whole boundary
        numSol(bNodes) = 0;
    case 2 % u=uex on whole boundary
        numSol(bNodes) = exSolNX(bNodes);
end

% ---------------------------------------------------------
% solution at nodes on the boundary
% ---------------------------------------------------------
solnbNodes = zeros(nNewNodes,1);
solebNodes = zeros(nNewNodes,1);
% solnbNodes(bNodes) = numSol(bNodes);
% solebNodes(bNodes) = exSolNX(bNodes);
solnbNodes(bNodeInt) = numSol(bNodeInt);
solebNodes(bNodeInt) = exSolNX(bNodeInt);
sol.nbNodes = solnbNodes;
sol.ebNodes = solebNodes;

%-----------------------------------------
% modification of F
F = F - A*numSol;
% LU factorization
numSol(iNodes) = A(iNodes,iNodes)\F(iNodes); % don't care nodes on boundary
% numSol(iNodes) = gmres(A(iNodes,iNodes),F(iNodes)); % GMRES factorization
% add to var
sol.num = numSol;



%% ========================================================
% NUMERICAL SOLUTION 
% in STANDARD FEM (just for plotting)
% VhSol_i = numSol_i for i in nodesOmg1NotOnGam or nodesOmg2NotAroundGam
% VhSol_i = numSol_k(i) for i in nodesOmg2CT
% VhSol_i = numSol_i+numSol_k(i) for i in nodesOnGam
% =========================================================
VhSol(nodesOmg1NotGamCTs) = numSol(nodesOmg1NotGamCTs);
VhSol(nodesOmg2NotInCTs) = numSol(nodesOmg2NotInCTs);
VhSol(nodesCTsInOmg2NotGam) = numSol(newNodes(nodesCTsInOmg2NotGam));
VhSol(nodesOnGamOfCTs) = numSol(newNodes(nodesOnGamOfCTs))...
                                + numSol(nodesOnGamOfCTs);
bNodesOnGam = intersect(nodesOnGam,bNodesStd); % nodes both on Gam and on boundary
VhSol(bNodesOnGam) = numSol(bNodesOnGam);
VhSol = VhSol'; % transform to column-array for pdesurf
% add to var
sol.Vh = VhSol;



%% ========================================================
% A PRIORI ERROR ESTIMATES
% error between numerical & exact solutions (column-array)
% =========================================================
errSol = exSolNX - numSol; 

% ---------------------------------------------------------
% L2
% ---------------------------------------------------------
mL2 = getMatrixL2(NCTs1,NCTs2,areaNCTs1,areaNCTs2,...
           iPs,CTs,typeCTs,nodeCTs,areaCTs,nodesCTsInOmg2OnGam,msh,pa);
errL2 = errSol' * mL2 * errSol;
errL2 = sqrt(errL2);
% add to var
err.L2=errL2;

% ---------------------------------------------------------
% H1 : We need to find ||k^{1/2}Grad e||_L2 first
% ---------------------------------------------------------
mL2grad = getMatrixL2Grad(NCTs1,NCTs2,nodesCTsInOmg2OnGam,...
                    CTs,areaChildCTs,areaCTs,msh,pa);
errL2G = errSol' * mL2grad * errSol;
errL2G = sqrt(errL2G);
% errH1 = errL2^2 + errL2G^2;
% errH1 = sqrt(errH1);
% % add to var
err.L2G=errL2G; 
% err.H1=errH1;

% ---------------------------------------------------------
% ENorm : given in Barrau's thesis (p.33)
% ---------------------------------------------------------
% lam*||[e]||
if ~isempty(CTs) % there exist cut-triangles
    mlNJU = getMatNormJumpU(CTs,iPs,msh,pa,cp.lambda);
    errlNJU = errSol' * mlNJU * errSol;
else % interface go through all nodes
    errlNJU = 0;
end
errlNJU = sqrt(errlNJU);
errENorm = errL2G^2 + errlNJU^2;
errENorm = sqrt(errENorm);
% add to var
err.ENorm=errENorm; 
 

%% ||[uh]|| % jump of numerical solution
ss=size(cp.lambda,2);
coef = 1+zeros(1,ss);
if ~isempty(CTs) % there exist cut-triangles
    mNJU = getMatNormJumpU(CTs,iPs,msh,pa,coef);
    errNJU = numSol' * mNJU * numSol;
else % interface go through all nodes
    errNJU = 0;
end
errNJU = sqrt(errNJU);
% add to var
err.NJU=errNJU;  

% ---------------------------------------------------------
% ||grad u-uex||_L2
% ---------------------------------------------------------
if ~isempty(CTs) % there exist cut-triangles
    mNJGU = getMatNormJumpGradU(CTs,iPs,msh,pa);
%     errNJGU = errSol' * mNJGU * errSol;
    errNJGU = numSol' * mNJGU * numSol;
else % interface go through all nodes
    errNJGU = 0;
end
errNJGU = sqrt(errNJGU);
% add to var
err.NJGU=errNJGU;


%% plot err |uex-uh|
% errSol in standard FE
errSolPlot = errSol; % uex-uh
% errSolPlot = abs(errSol); % |uex-uh|
errSolVh(nodesOmg1NotGamCTs) = errSolPlot(nodesOmg1NotGamCTs);
errSolVh(nodesOmg2NotInCTs) = errSolPlot(nodesOmg2NotInCTs);
errSolVh(nodesCTsInOmg2NotGam) = errSolPlot(newNodes(nodesCTsInOmg2NotGam));
errSolVh(nodesOnGamOfCTs) = errSolPlot(newNodes(nodesOnGamOfCTs))...
                                + errSolPlot(nodesOnGamOfCTs);
errSolVh(bNodesOnGam) = errSolPlot(bNodesOnGam);
errSolVh = errSolVh'; % transform to column-array for pdesurf

% add to var
sol.errSolVh=errSolVh;



%% ========================================================
% POSTERIORI ESTIMATES
% =========================================================

useit1=0;
if useit1
    % f in standard FEM 
    % ------------------------
    % (in barrau's model, f is the same for both subdomains but we will treat
    %   like considering the genreal case)
    defF = model.defF;
    fSD = zeros(nNodes,1); % column-array
    % fSD(nodesOmg1Gam) = defF(points(1,nodesOmg1Gam),...
    %         points(2,nodesOmg1Gam),1,pa);
    % fSD(nodesOmg2NotInCTs) = defF(points(1,nodesOmg2NotInCTs),...
    %         points(2,nodesOmg2NotInCTs),2,pa);
    % fSD(nodesCTsInOmg2NotGam) = defF(points(1,nodesCTsInOmg2NotGam),...
    %         points(2,nodesCTsInOmg2NotGam),2,pa);
    fSD(nodesOmg1) = defF(points(1,nodesOmg1),points(2,nodesOmg1),1,pa); % nodes inside Omg1
    fSD(nodesOnGam) = defF(points(1,nodesOnGam),points(2,nodesOnGam),1,pa); % nodes on Gam
    fSD(nodesOmg2) = defF(points(1,nodesOmg2),points(2,nodesOmg2),2,pa); % nodes inside Omg2

    % f in NXFEM 
    % ------------------------
    fNX = zeros(nNodes+nNewNodes,1); % column-array
    fNX(stdNodes) = fSD(stdNodes);
    fNX(newNodes(nodesInCTs)) = fSD(nodesInCTs);

    % etaK
    % -----------------------------
    metaK = getMatEstEtaK(NCTs1,NCTs2,areaNCTs1,areaNCTs2,...
                 iPs,CTs,typeCTs,nodeCTs,areaCTs,nodesCTsInOmg2OnGam,hT,pa,msh);
    erretaK = fNX'*metaK*fNX;
    erretaK = sqrt(erretaK);

    % etaS
    % -----------------------------
    if ~isempty(CTs) % there exist cut-triangles
        metaS = getMatEstEtaS(iPs,CTs,phi,pa,msh);
        erretaS = numSol'*metaS*numSol; 
    else % interface go through all nodes
        erretaS = 0;
    end
    erretaS = sqrt(erretaS);

    % zetaS
    %-----------------------------
    if ~isempty(CTs) % there exist cut-triangles
        lamZetaS = getLamZetaS(areaChildCTs,CTs,iPs,pa,msh);
        mzetaS = getMatNormJumpU(CTs,iPs,msh,pa,lamZetaS);
        errzetaS = numSol'*mzetaS*numSol;
    else % interface go through all nodes
        errzetaS = 0;
    end
    errzetaS = sqrt(errzetaS);

    % add to var
    %-----------------------------
    err.etaK=erretaK; err.etaS=erretaS; err.zetaS=errzetaS;
end % of useit1



%% ========================================================
% additional infor
% ========================================================
% plotting
pplot.gamh = getGamh(iPs);
pplot.iPs = iPs;
pplot.uNVCTs = uNVCTs;

end