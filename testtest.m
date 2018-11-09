% pdemesh(points,edges,triangles); % plot mesh
%     hold on;
    
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
       
    
    gvx = gvnew.x; gvy = gvnew.y;
    gvx(1,NCTs1(5,:))=0; gvy(1,NCTs1(5,:))=0;
    gvx(1,NCTs2(5,:))=0; gvy(1,NCTs2(5,:))=0;
    pdeplot(points,edges,triangles(1:3,:),'FlowData',[gvx; gvy]); 
    hold off;