%% This file is used to test ||phi0||_L2(Gam_h) and ||phi0-phi0_h||_L2(Omg)
% and compare with the results in freefem++ (cf. test_CR_phi0.edp)
% RESULT: - the same (with a very little different
%         - getNormL2fhf & getNormL2foGh are exact!!!!
%       ( for both regular and irregular meshes)


%% add path of functions
addpath(genpath('func')); % add all necessary functions

% nseg_array = [37, 57, 77, 101];
% nseg_array = 57;
nseg_array = zeros(1,5);

% hTmax = zeros(size(nseg_array,2),1);

hTmax = [0.1964185496	0.09701703355	0.04780764512	0.02450774691	0.01276380971];
hTmax = hTmax';

normPhionGh0 = zeros(size(nseg_array,2),1);
normPhih0PhionOmg = zeros(size(nseg_array,2),1);

order_phiGh = zeros(size(nseg_array,2),1);
order_phi_phih = zeros(size(nseg_array,2),1);


%% Fixed parameters
pa.degP1D = 3; % Gaussian quadrature points in 1D (polinomial functions)
pa.degP2D = 4; % Gaussian quadrature points in 2D (polinomial functions)
pa.degN = 8; % Gaussian quadrature points in 2D (non-polynomial functions)
% degree-#OfPoints : 1-1, 2quick -3, 3-4, 4-6, 5-7, 6-12, 7-13, 8-16, 9-19, 10-25, 11-27, 12-33
pa.tol = eps(1e3); % tolerance, 1e-14



%% models
% model = model_levelset_x; % Becker's test case with interface: x=x_0
model = model_levelset_vortex;  % Niklas' test case
GeoDom = model.domain(); % domain
velo = model.velo; % Velocity (grad of potential in other cases)
reguMesh = 1;

for iii=1:size(nseg_array,2)
% for iii=1:6

    %% setting
    nSeg = nseg_array(iii); % mesh settings

%     nSeg = 2^(iii+3);

    %% mesh setting up
%     if ~reguMesh % not regular mesh?
%         hEdgeMax = 2/nSeg;
%         [points,edges,triangles] = initmesh(GeoDom,'hmax',hEdgeMax); %irregular
%     else
%         [points,edges,triangles] = poimesh(GeoDom,nSeg,nSeg); % regular
%     end

%     file = strcat('Th_regular',num2str(iii-1),'.msh');
    file = strcat('Th_irregular',num2str(iii-1),'.msh');

    [points,edges,triangles] = getMeshFromFF(file);

    msh.p = points; msh.t = triangles; msh.e = edges; % save to msh
    x = points(1,:); % x-coordinate of points
    y = points(2,:); % y-coordinate of points
    % diameter (longest side) of each triangle: 1 x nTs
    msh.hT = getDiam(msh); % 1 x number of triangles
    msh.hTmax = max(msh.hT); % maximum of all diameters
%     hTmax(iii) = msh.hTmax;

    
    phi = model.defPhi(x,y,pa); % 1 x number of points (row array)
%     phi(abs(phi)<pa.tol) = 0; % find phi which are very small (~0) and set to 0

    tris = getTriangles(phi,msh,pa); % tris has 3 factors (structure var)
    CTs=tris.CTs;
    msh.nStd = size(points,2); % number of standard nodes
    CT = getInfoCTs(CTs,phi,msh,pa); % CT has many factors (structure var)
    iPs = CT.iPs;


    phi0 = phi;
    if ~isempty(CTs)
        normPhionGh0(iii) = getNormL2foGh(msh,pa,CTs,iPs,model.defPhi); % ||phi^0||_L2(Gam_h^0)
        normPhih0PhionOmg(iii) = getNormL2fhf(msh,pa,phi0,model.defPhi); % ||phi_h^0 - phi^0||_L2(Omg)
    end

end

for i=2:size(nseg_array,2)
    order_phiGh = log(normPhionGh0(i-1)/normPhionGh0(i))/log(hTmax(i-1)/hTmax(i));
    order_phi_phih = log(normPhih0PhionOmg(i-1)/normPhih0PhionOmg(i))/log(hTmax(i-1)/hTmax(i));
%     disp(hTmax(i));
    disp(order_phiGh);
    disp(order_phi_phih);
end