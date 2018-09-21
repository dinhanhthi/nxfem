function plotUNVwM(iPs,uNVCTs,msh)
% plot unit normal vector with mesh and gam_h
% Input: - mesh msh,
%        - unit normal vectors uNVCTs: 2 coor x nCTs 
%        - intersection points: 2 coor x 2 cut points x nCTs
%        - interface gamh: 
% Output: the plot

mp = getMidPointsCTs(iPs); % get the mid points (2 coor x nCTs)
nCTs = size(iPs,3); % number of cut triangles
% nf = size(get(0,'Children'),1);
% figure(nf+1);
figure(4);
% pdemesh(msh.p,msh.e,msh.t,'NodeLabels','off','ElementLabels','off'); % mesh
pdemesh(msh.p,msh.e,msh.t,'NodeLabels','on','ElementLabels','off'); % mesh
hold on
for t=1:nCTs
    plot(iPs(1,:,t),iPs(2,:,t),'-r','LineWidth',1); % gam_h on each cut triangle
end
quiver(mp(1,:),mp(2,:),uNVCTs(1,:),uNVCTs(2,:),0.5);
hold off
end