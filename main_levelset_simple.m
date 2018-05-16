%% =======================================================================
% This file is used to test with a very simple level set function
% ------------------------------------------------------------------------
% PURPOSE: Coding level set + NOT YET couple with NXFEM (only phi)
% After this work, test with Becker's model
% ------------------------------------------------------------------------


%% =======================================================================
% DOMAIN: [-1,1]X[-1,1]
%-------------------------------------------------------------------------
% MODELS: cf. model_barrau_x
%=========================================================================


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
model = model_barrau_x; % Becker's test case with interface: x=x_0
% typeG: modify directly in defG.m



%-------------------------------------------------------------------------
% Settings
%-------------------------------------------------------------------------
nSeg = 51; % mesh settings
velo = [.01;0]; % Velocity (grad of potential in other cases)
pa.xi = -0.21; % initial interface
pa.reguMesh = 0; % use regular mesh or not?
pa.smallCut = 0; % ignore small-support basis (1=ignore,0=no)
pa.tH = 100; % to find the small support using (20) or (21) in arnold 2008


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
% Level set function (INITIAL)
%-------------------------------------------------------------------------
phi = model.defPhi(x,y,pa); % 1 x number of points (row array)
phi(abs(phi)<pa.tol)=0; % find phi which are very small (~0) and set to 0




%% =======================================================================
% EACH TIME STEP
%=========================================================================
maxStep = 10; % max number of steps
dt = 0.02; % time step
for ns = 1:maxStep
    
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
    b3Nodes = bN.e3; % node on \pt\Omg_3 (top)
    
    
    
    
    %% ====================================================================
    % PLOT interface phi
    %======================================================================
    nf = 0; % reset every loop to be sure uh, vh plotted on the same figure
%     nf = plotNXFEM(msh,iPs,nf,'eleLabel','off','nodeLabel','off'); % only mesh
    plotNXFEM(msh,iPs,nf,phi,'withMesh',true,'title','phi','dim',2,'export',false); % phi
    
%     abc = waitforbuttonpress; % wait for click

%     nCTs = size(iPs,3);
%     for t=1:nCTs
%         plot(iPs(1,:,t),iPs(2,:,t),'-r','LineWidth',1);
%         hold on
%     end
%     hold off
    
    pause(0.01);
    
    %% ====================================================================
    % SOLVING phi (level set function)
    % standard finite element
    %======================================================================
    
    %----------------------------------------------------------------------
    % stiffness matrix for level set
    %----------------------------------------------------------------------
%     Aphi = getGMlsChopp07(msh,pa,vhStd,dt,0.5); % use only 1 file

    % use form like in Arnold Book p.221
%     del = getDells(msh,vhStd);
    del = max(velo);
    Eij = getMEls_test(msh,pa,velo,del);
    Hij = getMHls_test(msh,pa,velo,del,dt,0.5);
    
    Aphi = Eij + Hij;
    
    %----------------------------------------------------------------------
    % load vector for level set
    %----------------------------------------------------------------------
    AFphi = Eij - Hij;
    phi = phi'; % row to column
    Fphi = AFphi*phi;
    
    %----------------------------------------------------------------------
    % seek phi
    %----------------------------------------------------------------------
    phi = Aphi\Fphi; % update phi
    phi = phi'; % column to row
    
%     clear msh tris CTs NCTs1 NCTs2 nodeCTs iN bN bNodes b3Nodes cpU cpV
end % for ns







