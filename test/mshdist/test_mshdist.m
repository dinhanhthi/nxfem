%% MUST PUT THIS FILE IN /nxfem/
% Check if mshdist works well?



%% add path of functions
addpath(genpath('func')); % add all necessary functions
pa.degP1D = 3; % Gaussian quadrature points in 1D (polinomial functions)
pa.degP2D = 4; % Gaussian quadrature points in 2D (polinomial functions)
pa.degN = 8; % Gaussian quadrature points in 2D (non-polynomial functions)
pa.tol = eps(1e3); % tolerance, 1e-14



%% models
% model = model_levelset_x; % Becker's test case with interface: x=x_0
model = model_levelset_vortex;  % Niklas' test case
GeoDom = model.domain(); % domain
velo = model.velo; % Velocity (grad of potential in other cases)
nNormGP = 0;
nNormtxt = cell(1,6);
norm_gradphi = zeros(1,6);



%% setting
nSeg = 51; % mesh settings
reguMesh = 1; % use regular mesh or not?
useFMM = 1;
useSUPG = 1;
showPlot = 1; % show plot or not
pa.smallCut = 0;
pa.tH = 10;
t = 0;
dt = 0.6;



%% mesh setting up
if ~reguMesh % not regular mesh?
    hEdgeMax = 2/nSeg;
    [points,edges,triangles] = initmesh(GeoDom,'hmax',hEdgeMax); %irregular
else
    [points,edges,triangles] = poimesh(GeoDom,nSeg,nSeg); % regular
end
msh.p = points; msh.t = triangles; msh.e = edges; % save to msh
x = points(1,:); % x-coordinate of points
y = points(2,:); % y-coordinate of points
msh.hT = getDiam(msh); % 1 x number of triangles
msh.hTmax = max(msh.hT); % maximum of all diameters
hTmax = msh.hTmax;



%% Level set function (INITIAL)
phi = model.defPhi(x,y,pa); % 1 x number of points (row array)
phi(abs(phi)<pa.tol)=0; % find phi which are very small (~0) and set to 0



%% triangles info
tris = getTriangles(phi,msh,pa); % tris has 3 factors (structure var)
CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
msh.nStd = size(points,2); % number of standard nodes
CT = getInfoCTs(CTs,phi,msh,pa); % CT has many factors (structure var)
nodeCTs=CT.nodes; areaChildCTs=CT.areaChild;iPs=CT.iPs;



%% calc ||grad phi||
nNormGP = nNormGP+1;
norm_gradphi(nNormGP) = getNormL2Gstd(msh,phi);
nNormtxt{nNormGP} = '(initial)';



%% plot initial phi
if showPlot
    titlePlot1 = 'initial';
    nf = 0; % reset every loop to be sure uh, vh plotted on the same figure
    nf = plotNXFEM(msh,iPs,nf,phi,'withMesh',true,'title',titlePlot1,'dim',3,'export',false); % phi
end


% phi = phi.^3;
% if showPlot
%     titlePlot2 = 'modified';
%     nf = plotNXFEM(msh,iPs,nf,phi,'withMesh',true,'title',titlePlot2,'dim',3,'export',false); % phi
% end
% nNormGP = nNormGP+1;
% norm_gradphi(nNormGP) = getNormL2Gstd(msh,phi);
% nNormtxt{nNormGP} = '(phi^3)';



%% create command to run mshdist outside matlab
path_nxfem = '/home/thi/Dropbox/git/nxfem/'; % thi's local machine
path_phi = strcat(path_nxfem,'mshdist/');
call_mshdist = strcat({'mshdist'},{' '},{path_phi},'/phi_test'); % run in terminal
call_mshdist = cell2mat(call_mshdist);



% %% mshdist
% if useFMM
%     mshdist_w_mesh(msh,path_phi,'phi_test'); % export to .mesh
%     mshdist_w_sol(msh,phi,path_phi,'phi_test'); % export to .sol
%     system(call_mshdist); % run 'mshdist file/to/phi' (redistancing)
%     phi = mshdist_r_sol(phi,path_phi,'phi_test'); % update phi
%     
%     %% plot after mshdist
%     if showPlot
%         titlePlot3 = 'with mshdist';
%         plotNXFEM(msh,iPs,nf,phi,'withMesh',true,'title',titlePlot3,'dim',3,'export',false); % phi
%     end
%     
%     nNormGP = nNormGP+1;
%     norm_gradphi(nNormGP) = getNormL2Gstd(msh,phi);
%     nNormtxt{nNormGP} = '(phi^3 - FMM)';
% end



% %% print errors
% fprintf('norm_gradPhi (initial): %f,\n',norm_gradphi1);
% fprintf('norm_gradPhi (phi^3): %f,\n',norm_gradphi2);
% fprintf('norm_gradPhi (phi^3 - FMM): %f,\n',norm_gradphi3);



%% phi after solving
%========================================================================
t = t+dt;
vel_t = @(x,y,sub) velo(x,y,t,sub);



%% get triangles
tris = getTriangles(phi,msh,pa); % tris has 3 factors (structure var)
CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;



%% get cut triangles' info
CT = getInfoCTs(CTs,phi,msh,pa); % CT has many factors (structure var)
nodeCTs=CT.nodes; areaChildCTs=CT.areaChild;iPs=CT.iPs;



%% find small cut
if pa.smallCut
    [tris,CT] = findSmallPhi_after(msh,pa,phi,tris,CT);
    clear CTs NCTs NCTs2 nodeCTs areaChildCTs iPs; % just in case
    CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
    nodeCTs=CT.nodes; areaChildCTs=CT.areaChild;iPs=CT.iPs;
end



%% find nodes
msh.nNew = nodeCTs.n; % number of new nodes (nodes around the interface)
msh.ndof = msh.nNew + msh.nStd; % number of dofs
msh.newNodes = getNewNodes(nodeCTs.all,msh.nStd); % vector contaning new numbering of nodes around interface, column
msh.node = getNodes(tris,nodeCTs,msh,phi,pa); % get all nodes



%% boundary condition
[iN,bN] = getibNodes(msh);
bNodes = bN.all; iNodes=iN.all;
b3Nodes = bN.e3; % node on \pt\Omg_3 (top)



%% get del_T
if useSUPG
    del = getDellsT(msh,vel_t,1e-3,.5); % Arnold's book p.223
else
    del = zeros(1,size(msh.t,2)); % without SUPG
end
%     del = ones(1,size(msh.t,2))*hTmax/2; % just for testing    



%% stiffness matrix    
Eij = getMEls(msh,pa,vel_t,del);
Hij = getMHls(msh,pa,vel_t,del,dt,0.5);
Aphi = Eij + Hij;



%% load vector
AFphi = Eij - Hij;
phi = phi'; % row to column
Fphi = AFphi*phi;



%% get phi
phi = Aphi\Fphi; % update phi
phi = phi'; % column to row



%% ||gradPhi||
nNormGP = nNormGP+1;
norm_gradphi(nNormGP) = getNormL2Gstd(msh,phi);
nNormtxt{nNormGP} = '(with dt)';



%% plot after solving phi
if showPlot
    titlePlot4 = strcat('dt = ',num2str(dt),' (without mshdist)');
    nf = plotNXFEM(msh,iPs,nf,phi,'withMesh',true,'title',titlePlot4,'dim',3,'export',false); % phi
end



%% mshdist
if useFMM
    mshdist_w_mesh(msh,path_phi,'phi_test'); % export to .mesh
    mshdist_w_sol(msh,phi,path_phi,'phi_test'); % export to .sol
    system(call_mshdist); % run 'mshdist file/to/phi' (redistancing)
    phi = mshdist_r_sol(phi,path_phi,'phi_test'); % update phi
end



%% ||gradPhi||
nNormGP = nNormGP+1;
norm_gradphi(nNormGP) = getNormL2Gstd(msh,phi);
nNormtxt{nNormGP} = '(with dt - FMM)';



%% plot after solving phi + use FMM
if showPlot
    titlePlot4 = strcat('dt = ',num2str(dt),' (with mshdist)');
    plotNXFEM(msh,iPs,nf,phi,'withMesh',true,'title',titlePlot4,'dim',3,'export',false); % phi
end



%% print result
for i=1:nNormGP
   fprintf(strcat('norm_graphi ',nNormtxt{i},': %f,\n'),norm_gradphi(i)); 
end



%% get deltaT: 1 x nTs
function delT = getDellsT(msh,vel_t,eps,SD)
    % This function is only to be used type of velocity = (x,y,sub)
    % This function follows formula presented in Arnold's book, page 223 (and also cf. (7.14))
    % If in future, we need it more than once, I will put it in a separated file
    % Input: - msh.hT : h on each triangle : 1 x nTs
    %        - velo at each time step t: (x,y,sub)
    %        - msh.p to get all x,y coordinates
    %        - a control number SD = O(1) (cf. (7.14))
    %        - given small eps > 0 (cf. 7.14), at page 223, he took eps=1e-3
    % Output: a vector 1 x nTs

    x = msh.p(1,:); y = msh.p(2,:);

    % Take velo at each triangle
    veloPu = vel_t(x,y,1); % 1 x nPs
    veloPv = vel_t(x,y,2); % 1 x nPs
    veloTu = veloPu(msh.t(1:3,:)); % 3 x nTs
    veloTumax = max(abs(veloTu),[],1); % 1 x nTs
    veloTv = veloPv(msh.t(1:3,:)); % 3 x nTs
    veloTvmax = max(abs(veloTv),[],1); % 1 x nTs
    normVeloT = max(veloTumax,veloTvmax); % 1 x nTs

    normVeloT = max(eps,normVeloT); % 1 x nTs

    delT = SD*msh.hT./normVeloT;
end