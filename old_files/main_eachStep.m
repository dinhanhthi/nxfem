function [msh,err,sol,pplot,phi] = main_eachStep(hEdgeMax,GeoDom,pa,model)
% Find L2 err and H1 err for each step
% This funciton is used for finding convergence rate in main_CR.m
% Input: hEdgeMax, GeoDown to create a mesh
% Output: - error in L2 norm and H1 norm wrt this step
%         - some testing outputs
%         - some useful info about the mesh and solution



%% Get mesh info
if ~pa.reguMesh
    [points,edges,triangles] = initmesh(GeoDom,'hmax',hEdgeMax); % irregular
else
    numSeg = 2/hEdgeMax;
    [points,edges,triangles] = poimesh(GeoDom,numSeg,numSeg); % regular
end
msh.p = points; msh.t = triangles; msh.e = edges;
x = points(1,:); % x-coordinate of points
y = points(2,:); % y-coordinate of points
% diameter (longest side) of each triangle: 1 x nTs
msh.hT = getDiam(msh); % 1 x number of triangles
msh.hTmax = max(msh.hT); % maximum of all diameters



%% Level set function
phi = model.defPhi(x,y,pa); % 1 x number of points (row array)
phi(abs(phi)<pa.tol)=0; % find phi which are very small (~0) and set to 0



%% Get all triangles
tris = getTriangles(phi,msh,pa);
CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;



%% CT's info
CT = getInfoCTs(CTs,phi,msh,pa);
nodeCTs=CT.nodes; areaChildCTs=CT.areaChild; iPs=CT.iPs;
        


%% Small cut
if pa.smallCut
    [tris,CT] = findSmallPhi_after(msh,pa,phi,tris,CT);
    clear CTs NCTs NCTs2 nodeCTs areaChildCTs iPs; % just in case
    CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
    nodeCTs=CT.nodes; areaChildCTs=CT.areaChild;iPs=CT.iPs;
end



%% Nodes
msh.nNew = nodeCTs.n; % number of new nodes (nodes around the interface)
msh.nStd = size(points,2); % number of standard nodes
msh.ndof = msh.nNew + msh.nStd; % number of dofs
msh.newNodes = getNewNodes(nodeCTs.all,msh.nStd); % vector contaning new numbering of nodes around interface, column
msh.node = getNodes(tris,nodeCTs,msh,phi,pa); % get all nodes



%% Boundary nodes
[iN,bN] = getibNodes(msh);
bNodes = bN.all; iNodes=iN.all;



%% Exact solution in stdFEM
% exSol_i = exSol(x_i)
defExSol = model.defExSol;
uExStd = exInStd(defExSol,msh,pa);
sol.ex = uExStd;



%% Exact solution in NXFEM
% wExNX_i = wExStd_i for i is node of mesh
% wExNX_k(i) = wExStd_i for i in msh.node.CT.all
uExNX = interSTD2NX(uExStd,msh); % column array



%% Control paramaters
% depend on mesh
kap = model.kap(areaChildCTs,pa);
cp.kap1 = kap.kap1; cp.kap2 = kap.kap2; % kappa_i
cp.lambda = model.lam(pa,msh.hT,CTs); % penalty coef (not ghost penalty)
cp.kk1 = pa.kk1; cp.kk2 = pa.kk2;    % diff coef



%% Stiffness
% on all nodes including nodes on the boundary
A = getGMGG(tris,phi,CT,msh,pa,cp);



%% Load vector
defF = model.defF;
F = getLf(msh,pa,tris,CT,defF);


 
%% Applying boundary conditions
uhNX = zeros(msh.ndof,1); % column-array
typeBC = model.bc(); % get type of BCs
switch typeBC
    case 1 % u=o on whole boundary
        uhNX(bNodes) = 0;
    case 2 % u=uex on whole boundary
        uhNX(bNodes) = uExNX(bNodes);
end



%% Get solution
F = F - A*uhNX;
% LU factorization
uhNX(iNodes) = A(iNodes,iNodes)\F(iNodes); % don't care nodes on boundary
% uhNX(iNodes) = gmres(A(iNodes,iNodes),F(iNodes)); % GMRES factorization
sol.num = uhNX; % add to var



%% num sol in stdFEM
% VhSol_i = numSol_i for i in nodesOmg1NotOnGam or nodesOmg2NotAroundGam
% VhSol_i = numSol_k(i) for i in nodesOmg2CT
% VhSol_i = numSol_i+numSol_k(i) for i in msh.node.onG
uhStd = interNX2STD(uhNX,msh);
sol.Vh = uhStd; % add to var



%% Get errors
eU = uhNX - uExNX; % wex-wh
err.L2 = getNormL2(eU,tris,CT,msh,pa);

% cp1.kk1=1; cp1.kk2=1;
% cp1.kap1 = cp.kap1; cp1.kap2 = cp.kap2;
% err.L2G = getNormL2G(eU,tris,areaChildCTs,msh,cp1); % ||grad||_L2
err.L2G = getNormL2G(eU,tris,areaChildCTs,msh,cp); % ||kgrad||_L2
jumU = getNormJump1p2(eU,CTs,iPs,msh,pa);

% avegnU = getNormAveGn(eU,CTs,iPs,msh,cp1); % ||{gran w}||_{-1/2}
avegnU = getNormAveGn(eU,CTs,iPs,msh,cp); % ||{kgran w}||_{-1/2}
err.ENorm = err.L2G^2 + jumU^2+ avegnU^2;

% err.ENorm = err.L2^2 + err.L2G^2;

% coef = 4*pa.lamH*pa.kk1*pa.kk2/(pa.kk1+pa.kk2);
% err.ENorm = err.L2G^2 + coef*jumU^2;
% err.ENorm = err.L2G^2 + coef*jumU^2 + avegnU^2;

err.ENorm = sqrt(err.ENorm);



%% Errors (old)
% % =========================================================
% eU = uhNX - uExNX; % wex-wh
% err.L2 = getNormL2(eU,tris,CT,msh,pa);
% err.L2G = getNormL2G(eU,tris,areaChildCTs,msh,pa); % ||kgrad||_L2
% % cp1.kk1=1; cp1.kk2=1;
% % cp1.kap1 = cp.kap1; cp1.kap2 = cp.kap2;
% % err.L2G = getNormL2G(eU,tris,areaChildCTs,msh,cp1); % ||grad||_L2
% % ENorm : given in Barrau's thesis (p.33)
% % lam*||[e]||
% if ~isempty(CTs) % there exist cut-triangles
%     mlNJU = getMatNormJumpU(CTs,iPs,msh,pa,cp.lambda);
%     errlNJU = eU' * mlNJU * eU;
% else % interface go through all nodes
%     errlNJU = 0;
% end
% errlNJU = sqrt(errlNJU);
% errENorm = err.L2G^2 + errlNJU^2;
% errENorm = sqrt(errENorm);
% err.ENorm=errENorm; % add to var



%% just for plotting outside the box
pplot.gamh = getGamh(iPs);
pplot.iPs = iPs;
pplot.uNVCTs = CT.uN;


end