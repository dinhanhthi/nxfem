% (29-10-18) Last status is "wrong!!!" but don't know why, check
%   again in test_gradPhi.m
% Status: - with regular mesh: good
%         - irregular: good

clear all;

xDomVal = [0 1 1 0]; % x values of points constructing Omega
yDomVal = [0 0 1 1]; % corresponding y value
RectDom = [3,4,xDomVal,yDomVal]'; % rectangular domain "3" with "4" sides
GeoDom = decsg(RectDom);

reguMesh = 0;
nSeg = 15;
if ~reguMesh % not regular mesh?
    hEdgeMax = 2/nSeg;
    [points,edges,triangles] = initmesh(GeoDom,'hmax',hEdgeMax);    % irregular
else
    [points,edges,triangles] = poimesh(GeoDom,nSeg,nSeg);           % regular
end

% plot mesh
pdemesh(points,edges,triangles,'NodeLabels','on',...
                'ElementLabels','on');
            
msh.p = points; msh.t = triangles; msh.e = edges;   % save to msh

gradPhi = getGradPhi(msh.t,msh); % 2 coor x 3 vertices x nTris

% In the case irregular
% phi_i(x,y) = ax+by+c
i = 3; % between 1, 2, 3
t = 105;
switch i
    case 1
        j=2; k=3;
    case 2
        j=3; k=1;
    case 3
        j=1; k=2;
end
C = 1 - dot(gradPhi(:,i,t),msh.p(:,msh.t(i,t)));
phi_j = dot(gradPhi(:,i,t),msh.p(:,msh.t(j,t))) + C; % must be 0
phi_k = dot(gradPhi(:,i,t),msh.p(:,msh.t(k,t))) + C; % must be 0

disp(phi_j);
disp(phi_k);