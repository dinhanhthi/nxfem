% wanna know how pdegrad take at the center 
nSeg=2;
xDomVal = [0 1 1 0]; % x values of points constructing Omega
yDomVal = [0 0 1 1]; % corresponding y value
RectDom = [3,4,xDomVal,yDomVal]'; % rectangular domain "3" with "4" sides
GeoDom = decsg(RectDom);
[points,edges,triangles] = poimesh(GeoDom,nSeg,nSeg);

% plot mesh
pdemesh(points,edges,triangles,'NodeLabels','on',...
                'ElementLabels','on');
            
msh.p = points; msh.t = triangles; msh.e = edges;   % save to msh
u = ones(size(points,2),1); % sol=1

[vx,vy] = pdegrad(points, triangles, u);
% pdeplot(points, edges, triangles, 'FlowData', [vx', vy']);