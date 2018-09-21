function plotGradient(u,tris,msh,wM,nf)
% plot gradient of u on the triangles tris
% Input: - mesh msh
%        - u
%        - plot on which triangles? tris
%        - wM: plot with mesh or not? 1 or 0
%        - number of figure nf
% Output: figures

points=msh.p; triangles=msh.t; edges=msh.e;
ntris = size(tris,2); % number of triangles tris
gradU = getGradU(u,tris,msh); % get gradient of u on triangles tris
% 2 coor x 3 vertices x ntris


figure(nf)
for t=1:ntris
    xx = points(1,tris(1:3,t));
    yy = points(2,tris(1:3,t));
    gradux = gradU(1,:,t); % grad_x of u at 3 vertices of triangle t
    graduy = gradU(2,:,t); % grad_y of u at 3 vertices of triangle t
    quiver(xx,yy,gradux,graduy); % plot the velocity
    hold on
end

if wM
    figure(nf)
    pdemesh(points,edges,triangles,'NodeLabels','off','ElementLabels','off'); % plot mesh
else

hold off
end