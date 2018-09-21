function [msh,err,sol,pplot,phi] = main_eachStep(hEdgeMax,GeoDom,pa,model)
% Find L2 err and H1 err for each step
% IMPORTANT: this file contains etaK, etaS, zeta,....
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
msh.hT = getDiam(msh); % 1 x number of triangles
msh.hTmax = max(msh.hT); % maximum of all diameters
tol=pa.tol;
% level set function (interface)
phi = model.defPhi(x,y,pa); % 1 x number of points (row array)
phi(abs(phi)<tol)=0; % find values of phi which are very small (~0) and set them to 0
% exact solution function
defExSol = model.defExSol;



%% ========================================================
% GET INFORMATION OF TRIANGLES
% =========================================================

%-------------------------------------------------------------------------
% TRIANGLES
%-------------------------------------------------------------------------
[CTs,NCTs1,NCTs2] = getTriangles(phi,msh,pa);

%-------------------------------------------------------------------------
% ON CUT TRIANGLES
%-------------------------------------------------------------------------
[nodeCTs,typeCTs,areaChildCTs,iPs,uNVCTs] = getInfoCTs(CTs,phi,msh,pa);
        
%-------------------------------------------------------------------------
% Find small-cut triangles (idx in the OLD CTs)
%-------------------------------------------------------------------------
if pa.smallCut
    [tri2del,t2Omg1,t2Omg2] = findSmallPhi(msh,CTs,iPs,...
                                    nodeCTs,pa,typeCTs,phi);
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
        [nodeCTs,typeCTs,areaChildCTs,iPs,uNVCTs] = getInfoCTs(CTs,phi,msh,pa);
    end
end



%% ========================================================
% NODES
% =========================================================
msh.node.CT.onG = nodeCTs.onG; % nodes of CTs on Gam
msh.node.CT.iomg2 = nodeCTs.Omg2; % nodes in Omg2 (not on interface) of CTs
msh.node.CT.all = nodeCTs.all; % nodes in CTs, row-array
msh.nNew = nodeCTs.n; % number of new nodes (nodes around the interface)
msh.nStd = size(points,2); % number of standard nodes
msh.ndof = msh.nNew + msh.nStd; % number of dofs
msh.node.std = 1:msh.nStd; % standard nodes
msh.node.iomg1 = msh.node.std((phi<0)&(abs(phi)>tol)); 
    % all std nodes inside Omg1
msh.node.onG = msh.node.std(abs(phi)<tol); % all std nodes on Gam
msh.node.iomg2 = msh.node.std((phi>0)&(abs(phi)>tol)); 
    % all std nodes inside Omg2
msh.newNodes = getNewNodes(msh.node.CT.all,msh.nStd);
    % vector contaning new numbering of nodes around interface, column
msh.node.omg1 = unique(NCTs1(1:3,:)); 
msh.node.omg1 = union(msh.node.omg1,msh.node.CT.onG);
    % nodes inside Omg1 and on Gam, column-array
    % column-array, there still nodes on Gam in CTs which not belongs to 
    %   NCTs, that's why we need 2nd line
msh.node.omg2.all = unique(NCTs2(1:3,:)); 
    % column-array, there still nodes on Gam in CTs not belongs to NCTs
msh.node.omg2.all = union(msh.node.omg2.all,msh.node.CT.onG); % column-array
    % nodes inside Omg2 and on Gam
msh.node.omg2.notCT = setdiff(msh.node.omg2.all,msh.node.CT.all);
    % nodes in Omg2 but not in CTs (colum-array)
    % there may be still some nodes on the interface
msh.node.CT.omg2 = [msh.node.CT.iomg2,msh.node.CT.onG];
    % nodes in Omg2 and on G of CTs, row-array



%% ========================================================
% CONTROL PARAMETERS
% depend on mesh
% =========================================================
% cp = changing parameters
[cp.kap1,cp.kap2] = getKappa(areaChildCTs,pa); % follows Becker's guide
% cp.lambda = getLambda(areaChildCTs,CTs,iPs,pa,msh);
cp.lambda = getLambda(areaChildCTs,CTs,iPs,pa,msh);
cp.kk1 = pa.kk1; cp.kk2 = pa.kk2;


%% ========================================================
% STIFFNESS MATRIX
% on all nodes including nodes on the boundary
% =========================================================
A = getGMGG(NCTs1,NCTs2,CTs,phi,areaChildCTs,iPs,uNVCTs,msh,pa,cp);
% A in this case looks like Atw of tw equation (cf. main_sys_linda.m)



%% ========================================================
% LOAD VECTOR
% the right hand side
% =========================================================
defF = model.defF;
F = getLf(NCTs1,NCTs2,CTs,iPs,typeCTs,nodeCTs,msh,pa,defF);
% rhs in this case looks like rhs of w equation (cf. main_sys_linda.m)

 
%% ========================================================
% BOUNDARY CONDITIONS
% =========================================================
bNodesStd = unique([edges(1,:) edges(2,:)]); % standard boundary nodes,
% boundary nodes around interface
bNodesInCTs = intersect(bNodesStd,msh.node.CT.all);
bNodeInt = bNodesInCTs;
% NEW boundary nodes around interface
bNodesInCTs = msh.newNodes(bNodesInCTs); 
bNodeInt = [bNodeInt,bNodesInCTs']; % std and new nodes on boundary cap interface
bNodes = [bNodesStd,bNodesInCTs']; % all boundary nodes (std + new)
iNodesStd = setdiff(1:msh.nStd,bNodesStd); % standard inner nodes
% inner nodes around interface
iNodesInCTs = intersect(iNodesStd,msh.node.CT.all); 
% NEW boundary nodes around interface
iNodesInCTs = msh.newNodes(iNodesInCTs);
iNodes = [iNodesStd,iNodesInCTs']; % all inner nodes (std + new)



%% ========================================================
% EXACT SOLUTION
% in standard FEM space (just for plotting)
% exSol_i = exSol(x_i)
% =========================================================
exSol = zeros(msh.nStd,1); % column-array
exSol(msh.node.iomg1) = defExSol(points(1,msh.node.iomg1),points(2,msh.node.iomg1),1,pa); % nodes inside Omg1
exSol(msh.node.onG) = defExSol(points(1,msh.node.onG),points(2,msh.node.onG),1,pa); % nodes on Gam
exSol(msh.node.iomg2) = defExSol(points(1,msh.node.iomg2),points(2,msh.node.iomg2),2,pa); % nodes inside Omg2
% add to var
sol.ex = exSol;



%% ========================================================
% EXACT SOLUTION
% in NXFEM (for finding |exSol-sol|)
% exSolNX_i = exSol_i for i is node of mesh
% exSolNX_k(i) = exSol_i for i in nodesAroundGam
% =========================================================
exSolNX = zeros(msh.ndof,1); % column-array
exSolNX(msh.node.std) = exSol(msh.node.std);
exSolNX(msh.newNodes(msh.node.CT.all)) = exSol(msh.node.CT.all);



%% ========================================================
% PRECONDITIONER
% nodes outside CTs region
%% ========================================================
useit2=0;
if useit2 
    P = getMatPrecond(NCTs1,NCTs2,CTs,pa,areaChildCTs,msh);
end % useit2



%% ========================================================
% SEEK NUMERICAL SOLUTION
% in NXFEM space
% =========================================================
numSol = zeros(msh.ndof,1); % column-array
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
solnbNodes = zeros(msh.nNew,1);
solebNodes = zeros(msh.nNew,1);
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
% VhSol_i = numSol_i+numSol_k(i) for i in msh.node.onG
% =========================================================
VhSol(msh.node.omg1) = numSol(msh.node.omg1);
VhSol(msh.node.omg2.notCT) = numSol(msh.node.omg2.notCT);
VhSol(msh.node.CT.iomg2) = numSol(msh.newNodes(msh.node.CT.iomg2));
VhSol(msh.node.CT.onG) = numSol(msh.newNodes(msh.node.CT.onG))...
                                + numSol(msh.node.CT.onG);
bNodesOnGam = intersect(msh.node.onG,bNodesStd); % nodes both on Gam and on boundary
VhSol(bNodesOnGam) = numSol(bNodesOnGam);
VhSol = VhSol'; % transform to column-array for pdesurf
sol.Vh = VhSol; % add to var



%% ========================================================
% A PRIORI ERROR ESTIMATES
% error between numerical & exact solutions (column-array)
% =========================================================
errSol = exSolNX - numSol; 

% ---------------------------------------------------------
% L2
% ---------------------------------------------------------
mL2 = getMatrixL2(NCTs1,NCTs2,iPs,CTs,typeCTs,nodeCTs,msh,pa);
errL2 = errSol' * mL2 * errSol;
errL2 = sqrt(errL2);
% add to var
err.L2=errL2;



% ---------------------------------------------------------
% H1 : We need to find ||k^{1/2}Grad e||_L2 first
% ---------------------------------------------------------
err.L2 = getErrL2(errSol,NCTs1,NCTs2,CTs,iPs,typeCTs,nodeCTs,msh,pa);

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
errSolVh(msh.node.omg1) = errSolPlot(msh.node.omg1);
errSolVh(msh.node.omg2.notCT) = errSolPlot(msh.node.omg2.notCT);
errSolVh(msh.node.CT.iomg2) = errSolPlot(msh.newNodes(msh.node.CT.iomg2));
errSolVh(msh.node.CT.onG) = errSolPlot(msh.newNodes(msh.node.CT.onG))...
                                + errSolPlot(msh.node.CT.onG);
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
    fSD = zeros(msh.nStd,1); % column-array
    % fSD(msh.node.omg1) = defF(points(1,msh.node.omg1),...
    %         points(2,msh.node.omg1),1,pa);
    % fSD(msh.node.omg2.notCT) = defF(points(1,msh.node.omg2.notCT),...
    %         points(2,msh.node.omg2.notCT),2,pa);
    % fSD(msh.node.CT.iomg2) = defF(points(1,msh.node.CT.iomg2),...
    %         points(2,msh.node.CT.iomg2),2,pa);
    fSD(msh.node.iomg1) = defF(points(1,msh.node.iomg1),points(2,msh.node.iomg1),1,pa); % nodes inside Omg1
    fSD(msh.node.onG) = defF(points(1,msh.node.onG),points(2,msh.node.onG),1,pa); % nodes on Gam
    fSD(msh.node.iomg2) = defF(points(1,msh.node.iomg2),points(2,msh.node.iomg2),2,pa); % nodes inside Omg2

    % f in NXFEM 
    % ------------------------
    fNX = zeros(msh.ndof,1); % column-array
    fNX(msh.node.std) = fSD(msh.node.std);
    fNX(newNodes(msh.node.CT.all)) = fSD(msh.node.CT.all);

    % etaK
    % -----------------------------
    metaK = getMatEstEtaK(NCTs1,NCTs2,iPs,CTs,typeCTs,nodeCTs,pa,msh);
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