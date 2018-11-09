function plotInterface(msh,pa,phi,iPs)
% Just for plotting the interface
% For the more complicated plot, cf. plotNXFEM.m


points=msh.p;

% plot intersections
nCTs = size(iPs,3);
for it=1:nCTs
    plot(iPs(1,:,it),iPs(2,:,it),'-r','LineWidth',1);
    hold on
end

% plot segments Gh conciding to mesh's edges
segment = getMeshSegmentOnGh(msh,pa,phi);
for it=1:size(segment,2)
    plot(points(1,segment(:,it)),points(2,segment(:,it)),'-r','LineWidth',1);
    hold on
end


end