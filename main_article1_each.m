function [hTmax,errW,errV] = main_article1_each(nSeg)
%% This files used in each step for finding convergence rate
% file to run: main_article1.m


%% =======================================================================
% PARAMETERS
% Note that parameters are different for equations of w and u
%=========================================================================

%-------------------------------------------------------------------------
% Fixed parameters
%-------------------------------------------------------------------------
pa.degP1D = 3; % Gaussian quadrature points in 1D (polinomial functions)
pa.degP2D = 4; % Gaussian quadrature points in 2D (polinomial functions)
pa.degN = 8; % Gaussian quadrature points in 2D (non-polynomial functions)
% degree-#OfPoints : 1-1, 2-3, 3-4, 4-6, 5-7, 6-12,
%                    7-13, 8-16, 9-19, 10-25, 11-27, 12-33
pa.tol = eps(1e3); % tolerance, 1e-14
model = model_article1; % choose model. cf. file model_article1.m
% typeG: modify directly in defG.m

%-------------------------------------------------------------------------
% General settings
%-------------------------------------------------------------------------
wannaPlot = 0; % wanna plot?
useNewton = 0; % wanna use Newton method for finding v or not?

%-------------------------------------------------------------------------
% Deleting small cut
%-------------------------------------------------------------------------
pa.smallCut = 0; % ignore small-support basis (1=ignore,0=no)
pa.tH = 1e2; % to find the small support using (20) or (21) in arnold 2008

%-------------------------------------------------------------------------
% Penalty parameters (\gam\int [w][varphi])
%-------------------------------------------------------------------------
pa.lamHw = 1e11; % penalty coefficient for w
pa.lamHv = 1e8; % penalty coefficient for v

%-------------------------------------------------------------------------
% Ghost penalty
%-------------------------------------------------------------------------
pa.useGP = 0; % wanna use ghost penalty term?
pa.gam1 = 1e-7; % parameter for 1st term
pa.gam2 = 1e-7 ; % parameter for 2nd term

%-------------------------------------------------------------------------
% Mesh settings
%-------------------------------------------------------------------------
pa.reguMesh = 0; % use regular mesh or not?

%-------------------------------------------------------------------------
% Model parameters
%-------------------------------------------------------------------------
pa.alp1 = 1; pa.alp2 = 100; % diff coeff for eqn w
pa.bet1 = 0.5; pa.bet2 = 100; % diff coeff for eqn v
% (bet2 doesn't take affect because v=0 in Omg2)
pa.r0 = 0.6; % interface
pa.lamSys = 1; % coef lam in system settings




%% =======================================================================
% DOMAIN
%=========================================================================
GeoDom = model.domain(); % domain

%-------------------------------------------------------------------------
% Mesh setting up
%-------------------------------------------------------------------------
if ~pa.reguMesh % not regular mesh?
    hEdgeMax = 2/nSeg;
    [points,edges,triangles] = initmesh(GeoDom,'hmax',hEdgeMax); %irregular
else
    [points,edges,triangles] = poimesh(GeoDom,nSeg,nSeg); % regular
end
msh.p = points; msh.t = triangles; msh.e = edges; % save to msh
x = points(1,:); % x-coordinate of points
y = points(2,:); % y-coordinate of points
% diameter (longest side) of each triangle: 1 x nTs
msh.hT = getDiam(msh); % 1 x number of triangles
msh.hTmax = max(msh.hT); % maximum of all diameters

hTmax = msh.hTmax;

%-------------------------------------------------------------------------
% Level set function
%-------------------------------------------------------------------------
phi = model.defPhi(x,y,pa); % 1 x number of points (row array)
phi(abs(phi)<pa.tol)=0; % find phi which are very small (~0) and set to 0



%% =======================================================================
% GET INFORMATION OF TRIANGLES
% The same for both equations of tw and u
%=========================================================================

%-------------------------------------------------------------------------
% Triangles
%-------------------------------------------------------------------------
tris = getTriangles(phi,msh,pa); % tris has 3 factors (structure var)
CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;

%-------------------------------------------------------------------------
% On cut triangles
%-------------------------------------------------------------------------
CT = getInfoCTs(CTs,phi,msh,pa); % CT has many factors (structure var)
nodeCTs=CT.nodes; areaChildCTs=CT.areaChild;iPs=CT.iPs;
        
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
        clear nodeCTs typeCTs areaChildCTs iPs uNVCTs areaCTs; % in case
        CT = getInfoCTs(CTs,phi,msh,pa);
        nodeCTs=CT.nodes; areaChildCTs=CT.areaChild;iPs=CT.iPs;
    end
end



%% =======================================================================
% NODES
%=========================================================================
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

% w
%-------------------------------------------------------------------------
defExSol = model.defWex;
wExStd = exInStd(defExSol,msh,pa);

% u
%-------------------------------------------------------------------------
defExSol = model.defUex;
uExStd = exInStd(defExSol,msh,pa);

% v
%-------------------------------------------------------------------------
defExSol = model.defVex;
vExStd = exInStd(defExSol,msh,pa);



%% ========================================================
% EXACT SOLUTION in NXFEM
% wExNX_i = wExStd_i for i is node of mesh
% wExNX_k(i) = wExStd_i for i in msh.node.CT.all
% =========================================================

% w
%-------------------------------------------------------------------------
wExNX = interSTD2NX(wExStd,msh); % column array

% v
%-------------------------------------------------------------------------
vExNX = interSTD2NX(vExStd,msh); % column array

% u
%-------------------------------------------------------------------------
uExNX = interSTD2NX(uExStd,msh); % column array



%% =======================================================================
% CONTROL PARAMETERS
% depend on mesh and different for w and u
% in child-functions, it's the variable "cp"
%=========================================================================
kapW = model.kapW(areaChildCTs,pa);
cpW.kap1 = kapW.kap1; cpW.kap2 = kapW.kap2; % kappa_i
cpW.kk1 = pa.alp1; cpW.kk2 = pa.alp2;    % diff coef
cpW.lambda = model.lamW(cpW,msh.hT,CTs,pa); % penalty coef (not ghost penalty)

kapV = model.kapV(areaChildCTs,pa);
cpV.kap1 = kapV.kap1; cpV.kap2 = kapV.kap2; % kappa_i
cpV.kk1 = pa.bet1; cpV.kk2 = pa.bet2;    % diff coef
cpV.lambda = model.lamV(cpV,msh.hT,CTs,pa); % penalty coef (not ghost penalty)




%% =======================================================================
% SOLVING W
%=========================================================================

%-------------------------------------------------------------------------
% Stiffness matrix (all nodes including nodes on boundary)
%-------------------------------------------------------------------------
Aw = getGMGG(tris,phi,CT,msh,pa,cpW);

%-------------------------------------------------------------------------
% Load vector (all nodes including nodes on boundary)
%-------------------------------------------------------------------------
defFw = model.defFw;
Fw = getLf(tris,CT,msh,pa,defFw);

%-------------------------------------------------------------------------
% BCs
%-------------------------------------------------------------------------
whNX = zeros(msh.ndof,1); % column-array
typeBC = model.bcW(); % get type of BCs
switch typeBC
    case 1 % u=o on whole boundary
        whNX(bNodes) = 0;
    case 2 % u=uex on whole boundary
        whNX(bNodes) = wExNX(bNodes);
end

%-------------------------------------------------------------------------
% Solving w
%-------------------------------------------------------------------------
Fw = Fw - Aw*whNX; % modification of F
% LU factorization
whNX(iNodes) = Aw(iNodes,iNodes)\Fw(iNodes); % don't care nodes on boundary
% whNX(iNodes) = gmres(Aw(iNodes,iNodes),Fw(iNodes)); % GMRES factorization



%% =======================================================================
% SOLVING V
%=========================================================================

% initial solution
% ------------------------------------------------------------------------
vold = zeros(msh.ndof,1);
% vold = ones(msh.ndof,1);
% vold = vExNX-0.05;

% tolerance
% ------------------------------------------------------------------------
itol = 1e-6;
delL2 = getNormL2(vold - vExNX,tris,CT,msh,pa);
Vip1L2 = getNormL2(vExNX,tris,CT,msh,pa); % |v_{i+1}|_L2
difv = delL2/Vip1L2; % |del|_L2/|v_i+1|_L2
% fprintf('difv: %0.7f\n',difv);
% fprintf('difv: %0.7f\n',difv);

imax = 50; % maximum number of steps
step = 0;
defFv = model.defFv;
typeBC = model.bcV(); % get type of BCs
wS = getUold(whNX,msh); % w in each subdomain
while (difv > itol) && (step<imax)
    step = step+1;
    
    % analyze numSolu into each subdomain
    % --------------------------------------------------------------------
    voldEach = getUold(vold,msh);
    
    if ~useNewton % don't wanna use Newton method
        
        % global matrix
        % -----------------------------------------------------------------
        Av = getGMvAA(tris,phi,voldEach,wS,CT,msh,pa,cpV,cpW);
        
        % load vector
        % -----------------------------------------------------------------
        Fv = getLf(tris,CT,msh,pa,defFv);
        
        % bc
        % -----------------------------------------------------------------
        vnew = zeros(msh.ndof,1); % zero initial uh for each step
        switch typeBC
            case 1 % u=o on whole boundary
                vnew(bNodes) = 0;
            case 2 % u=uex on whole boundary
                vnew(bNodes) = vExNX(bNodes);
        end
        
        % solving for u
        % -----------------------------------------------------------------
        Fv = Fv - Av*vnew; % modification of F
        % LU factorization
        vnew(iNodes) = Av(iNodes,iNodes)\Fv(iNodes); % don't care nodes on boundary
        % uh(iNodes) = gmres(Au(iNodes,iNodes),Fu(iNodes)); % GMRES factorization
        
        del = vnew - vold;
        vold = vnew; % update for the next step
        
    else % use Newton method (DF(u)*del = F(u), solve for del)
        
        % DF(u)*del
        % -----------------------------------------------------------------
        Adel = getGMvAANewton(tris,phi,voldEach,wS,CT,msh,pa,cpV,cpW);
               
        % F(u)
        % -----------------------------------------------------------------
        
        Av = getGMvAA(tris,phi,voldEach,wS,CT,msh,pa,cpV,cpW); 
                    % like Av in normal iterative method
        Fv = getLf(tris,CT,msh,pa,defFv);
                    % like Fu in normal iterative method
        Fdel = Av*vold - Fv;
        
        % bc for del
        % -----------------------------------------------------------------
        del = zeros(msh.ndof,1); % column-array
        del(bNodes) = 0; % always
        
        % solve for del
        % -----------------------------------------------------------------
        Fdel = Fdel - Adel*del; 
        % LU factorization
        del(iNodes) = Adel(iNodes,iNodes)\Fdel(iNodes); % don't care nodes on boundary
        % del(iNodes) = gmres(Adel(iNodes,iNodes),Fdek(iNodes)); % GMRES factorization
        
        vold(bNodes) = vExNX(bNodes); % u=uex on bc
        vold(iNodes) = vold(iNodes) - del(iNodes); % u_i+1 = u_i - del, update for the next step
    end
    
    delL2 = getNormL2(del,tris,CT,msh,pa);
    Vip1L2 = getNormL2(vold,tris,CT,msh,pa);
    difv = delL2/Vip1L2; % |del|_L2/|u_i+1|_L2
%     difv = delL2;
%     fprintf('difv: %0.18f\n',difv);
end
vhNX = vold;



%% =======================================================================
% SOLVING U
% u = w - bet/(alp*lamSys)*v
% Note that v=0 on nodes on the interface
%=========================================================================
uhNX = getUh(whNX,vhNX,pa.bet1/(pa.lamSys*pa.alp1),...
                    pa.bet2/(pa.lamSys*pa.alp2),msh);



%% =======================================================================
% ERRORS
%=========================================================================
% ---------------------------------------------
% for w, |||_1
% ---------------------------------------------
eW = whNX - wExNX; % wex-wh
errW.L2 = getNormL2(eW,tris,CT,msh,pa);

cp1w.kk1=1; cp1w.kk2=1;
cp1w.kap1 = cpW.kap1; cp1w.kap2 = cpW.kap2;
errW.L2G = getNormL2G(eW,tris,areaChildCTs,msh,cp1w); % ||grad||_L2
% errW.L2G = getNormL2G(eW,tris,areaChildCTs,msh,cpW); % ||kgrad||_L2
jumW = getNormJump1p2(eW,CTs,iPs,msh,pa); % |[u]|_1/2

avegnW = getNormAveGn(eW,CTs,iPs,msh,cp1w); % |{grad v}|_{-1/2}
errW.ENorm = errW.L2G^2 + jumW^2+ avegnW^2;

% coef = 4*pa.lamHw*cpW.kk1*cpW.kk2/(cpW.kk1+cpW.kk2);
% errW.ENorm = errW.L2G^2 + coef*jumW^2;

errW.ENorm = sqrt(errW.ENorm);


% ---------------------------------------------
% for v, |||_2
% ---------------------------------------------
eV = vhNX - vExNX; % vex-vh
errV.L2 = getNormL2(eV,tris,CT,msh,pa);

cp1v.kk1=1; cp1v.kk2=1;
cp1v.kap1 = cpV.kap1; cp1v.kap2 = cpV.kap2;
errV.L2G = getNormL2G(eV,tris,areaChildCTs,msh,cp1v);  % ||grad||_L2
% errV.L2G = getNormL2G(eW,tris,areaChildCTs,msh,cpV); % ||kgrad||_L2
jumV = getNormJump1p2(eV,CTs,iPs,msh,pa); % |[u]|_1/2
aveV = getNormAveV(eV,CTs,iPs,msh,cpV,pa); % |{v}|_L2(Gam)

% option 1
avegnV = getNormAveGn(eV,CTs,iPs,msh,cp1v); % |{grad v}|_{-1/2}
% avegnV = getNormAveGn(eV,CTs,iPs,msh,cpV); % |{kgrad v}|_{-1/2}
errV.ENorm = errV.L2G^2 + jumV^2 + avegnV^2 + aveV^2;

% % option 2
% coef = 4*pa.lamHv*cpV.kk1*cpV.kk2/(cpV.kk1+cpV.kk2);
% errV.ENorm = errV.L2G^2 + coef*jumV^2 + aveV^2;

errV.ENorm = sqrt(errV.ENorm);



%% ========================================================
% NUMERICAL SOLUTION 
% in STANDARD FEM (just for plotting)
% whStd_i = whNX_i for i in msh.node.omg1 or msh.node.omg2.notCT
% whStd_i = whNX_k(i) for i in msh.node.CT.iomg2
% whStd_i = whNX_i+whNX_k(i) for i in msh.node.CT.onG
% =========================================================

% w
%-------------------------------------------------------------------------
whStd = interNX2STD(whNX,msh);

% v
%-------------------------------------------------------------------------
vhStd = interNX2STD(vhNX,msh);

% u
%-------------------------------------------------------------------------
uhStd = interNX2STD(uhNX,msh);


%% =======================================================================
% PLOTTING
%=========================================================================

if wannaPlot
    nf = 0; % used for plot new figure's window
%     nf = plotNXFEM(msh,iPs,nf,'eleLabel','off','nodeLabel','on'); % only mesh
    nf = plotNXFEM(msh,iPs,nf,wExStd,'withMesh',true,'title','wex',...
                'dim',3,'export',false); % wex
    nf = plotNXFEM(msh,iPs,nf,whStd,'withMesh',true,...
                    'title','wh','dim',3,'export',false); % wh
    nf = plotNXFEM(msh,iPs,nf,vExStd,'withMesh',true,'title','vex',...
                'dim',3,'export',false); % vex
    nf = plotNXFEM(msh,iPs,nf,vhStd,'withMesh',true,...
                    'title','vh','dim',3,'export',false); % vh
%     nf = plotNXFEM(msh,iPs,nf,uExStd,'withMesh',true,'title','uex',...
%                     'dim',3,'export',false); % uex
%     nf = plotNXFEM(msh,iPs,nf,uhStd,'withMesh',true,...
%                     'title','uh','dim',3); % uh

end



%% =======================================================================
% Variables
% ------------------------------------------------------------------------
% wExStd,uExStd,vExStd: exact solution of w,u,v in standard FEM space
% wExNX,uExNX,vExNX: exact solution of w,u,v in NXFEM
% whStd,uhStd,vhStd: numerical solution of w,u,v in standard FEM space
% whNX,uhNX,vhNX: numerical solition of w,u,v in NXFEM
% ------------------------------------------------------------------------
% errW.L2, errV.L2: |wex-wh|_L2
%=========================================================================

end
