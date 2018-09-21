pdemesh(points,edges,triangles,'NodeLabels','on','ElementLabels','on');
nCTs = size(iPs,3);
hold on;
for t=1:nCTs
    plot(iPs(1,:,t),iPs(2,:,t),'-r','LineWidth',1); % gam_h on each cut triangle
end
hold off;
