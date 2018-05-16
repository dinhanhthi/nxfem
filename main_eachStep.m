
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
    [points,edges,triangles] = poimesh(GeoDom,numSeg,numSeg); % regular
end
% save to msh
msh.p = points; msh.t = triangles; msh.e = edges;
x = points(1,:); % x-coordinate of points
y = points(2,:); % y-coordinate of points
% diameter (longest side) of each triangle: 1 x nTs
msh.hT = getDiam(msh); % 1 x number of triangles
msh.hTmax = max(msh.hT); % maximum of all diameters



%-------------------------------------------------------------------------
% Level set function
%-------------------------------------------------------------------------
phi = model.defPhi(x,y,pa); % 1 x number of points (row array)
phi(abs(phi)<pa.tol)=0; % find phi which are very small (~0) and set to 0



%% ========================================================
% GET INFORMATION OF TRIANGLES
% =========================================================

%-------------------------------------------------------------------------
% TRIANGLES
%-------------------------------------------------------------------------
tris = getTriangles(phi,msh,pa);
CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;

%-------------------------------------------------------------------------
% ON CUT TRIANGLES
%-------------------------------------------------------------------------
CT = getInfoCTs(CTs,phi,msh,pa);
nodeCTs=CT.nodes; areaChildCTs=CT.areaChild; iPs=CT.iPs;
        
%-------------------------------------------------------------------------
% Find small-cut triangles (idx in the OLD CTs)
%-------------------------------------------------------------------------
if pa.smallCut
    [tri2del,t2Omg1,t2Omg2] = findSmallPhi(msh,CTs,CT,pa,phi);
    % If there are small-cut triangles, remove them from CTs!!
    if ~isempty(tri2del)
        nCTs = size(CTs,2); % number of OLD cut triangles
        % get NEW not-cut triangles
        if ~isempty(t2Omg1)
            NCTs1 = [NCTs1,CTs(:,t2Omg1)]; % add more triangles to NCTs1
            tris.NCTs1=NCTs1;
        end
        if ~isempty(t2Omg2)
            NCTs2 = [NCTs2,CTs(:,t2Omg2)]; % add more triangles to NCTs2
            tris.NCTs2=NCTs2;
        end
        % get NEW cut triangles
        CTs = CTs(:,setdiff(1:nCTs,tri2del));
        tris.CTs=CTs;
        % find again all information
        clear nodeCTs typeCTs areaChildCTs iPs uNVCTs; % in case
        CT = getInfoCTs(CTs,phi,msh,pa);
        nodeCTs=CT.nodes; areaChildCTs=CT.areaChild; iPs=CT.iPs;
    end
end



%% ========================================================
% NODES
% =========================================================
msh.nNew = nodeCTs.n; % number of new nodes (nodes around the interface)
msh.nStd = size(points,2); % number of standard nodes
msh.ndof = msh.nNew + msh.nStd; % number of dofs
msh.newNodes = getNewNodes(nodeCTs.all,msh.nStd);
    % vector contaning new numbering of nodes around interface, column
msh.node = getNodes(tris,nodeCTs,msh,phi,pa); % get all nodes


%-------------------------------------------------------------------------
% boundary nodes and inner nodes
%-------------------------------------------------------------------------
[iN,bN] = getibNodes(msh);
bNodes = bN.all; iNodes=iN.all;


%% ========================================================
% EXACT SOLUTION
% in standard FEM space (just for plotting)
% exSol_i = exSol(x_i)
% =========================================================
defExSol = model.defExSol;
uExStd = exInStd(defExSol,msh,pa);
sol.ex = uExStd;


%% ========================================================
% EXACT SOLUTION in NXFEM
% wExNX_i = wExStd_i for i is node of mesh
% wExNX_k(i) = wExStd_i for i in msh.node.CT.all
% =========================================================
uExNX = interSTD2NX(uExStd,msh); % column array




%% ========================================================
% CONTROL PARAMETERS
% depend on mesh
% =========================================================
kap = model.kap(areaChildCTs,pa);
cp.kap1 = kap.kap1; cp.kap2 = kap.kap2; % kappa_i
cp.lambda = model.lam(pa,msh.hT,CTs); % penalty coef (not ghost penalty)
cp.kk1 = pa.kk1; cp.kk2 = pa.kk2;    % diff coef



%% ========================================================
% STIFFNESS MATRIX
% on all nodes including nodes on the boundary
% =========================================================
A = getGMGG(tris,phi,CT,msh,pa,cp);
% A in this case looks like Atw of tw equation (cf. main_sys_linda.m)



%% ========================================================
% LOAD VECTOR
% the right hand side
% =========================================================
defF = model.defF;
F = getLf(tris,CT,msh,pa,defF);
% rhs in this case looks like rhs of w equation (cf. main_sys_linda.m)

 
%% ========================================================
% BOUNDARY CONDITIONS
% =========================================================
uhNX = zeros(msh.ndof,1); % column-array
typeBC = model.bc(); % get type of BCs
switch typeBC
    case 1 % u=o on whole boundary
        uhNX(bNodes) = 0;
    case 2 % u=uex on whole boundary
        uhNX(bNodes) = uExNX(bNodes);
end


%% ========================================================
% SOLVING
% =========================================================
F = F - A*uhNX;
% LU factorization
uhNX(iNodes) = A(iNodes,iNodes)\F(iNodes); % don't care nodes on boundary
% uhNX(iNodes) = gmres(A(iNodes,iNodes),F(iNodes)); % GMRES factorization
sol.num = uhNX; % add to var



%% ========================================================
% NUMERICAL SOLUTION 
% in STANDARD FEM (just for plotting)
% VhSol_i = numSol_i for i in nodesOmg1NotOnGam or nodesOmg2NotAroundGam
% VhSol_i = numSol_k(i) for i in nodesOmg2CT
% VhSol_i = numSol_i+numSol_k(i) for i in msh.node.onG
% =========================================================
uhStd = interNX2STD(uhNX,msh);
sol.Vh = uhStd; % add to var


%% =======================================================================
% ERRORS
%=========================================================================
eU = uhNX - uExNX; % wex-wh
err.L2 = getNormL2(eU,tris,CT,msh,pa);

cp1.kk1=1; cp1.kk2=1;
cp1.kap1 = cp.kap1; cp1.kap2 = cp.kap2;
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



% %% ========================================================
% % ERRORS (old)
% % ENorm looks like in Barrau thesis
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



%% ========================================================
% additional infor
% ========================================================
% plotting
pplot.gamh = getGamh(iPs);
pplot.iPs = iPs;
pplot.uNVCTs = CT.uN;

end