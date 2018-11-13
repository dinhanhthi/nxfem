function plotNXFEM(nS,eS,mshs,lsf,msh,iPs)
% Containing all plots
% Input: - numerical solution nS & exact solution eS in standard FEM
%        - options for plotting mesh: mshs
%        - options for plotting level set function: lsf
% Output: figures depending on what users wanna see

nfg = 1; % number of figure
triangles=msh.t; points=msh.p; edges=msh.e;

%% plot mesh
% ---------------------------------------------------------------
    function plotMesh(yOn,msh) % plot solution with mesh or not?
       if yOn
           hold on
           pdemesh(msh.p,msh.e,msh.t); % plot mesh
           hold off 
       end
    end
% ---------------------------------------------------------------

%% plot interface gamma_h
    function plotGamh(iPs)
        hold on
        nCTs = size(iPs,3);
        for t=1:nCTs
            plot(iPs(1,:,t),iPs(2,:,t),'-r','LineWidth',1); % gam_h on each cut triangle
        end
        hold off
    end


%% ========================================================
% MESH
% =========================================================
if mshs{1}
    figure(nfg); nfg = nfg+1;
    pdemesh(points,edges,triangles,'NodeLabels',mshs{2},'ElementLabels',mshs{3});
    plotGamh(iPs); % plot gam_h
end


%% ========================================================
% LEVEL SET FUNCTION
% =========================================================
if lsf{2} % user wanna plot level set function
    figure(nfg); nfg = nfg+1;
    if lsf{3}==2 % 2D
        pdeplot(points,edges,triangles,'XYData',lsf{1},'Title','Level Set function');
        plotMesh(lsf{4},msh); % plot with mesh
        plotGamh(iPs); % plot gam_h
    elseif lsf{3}==3 % 3D
        pdeplot(points,[],triangles,'XYData',lsf{1},'XYStyle','interp',...
         'ZData',lsf{1},'ZStyle','continuous',...
         'ColorBar','off','Title','Level set func');  
        plotMesh(lsf{4},msh); % plot with mesh
    else % both 1D and 2D
        pdeplot(points,edges,triangles,'XYData',lsf{1},'Title','Numerical solution');
        plotMesh(lsf{4},msh); % plot with mesh
        plotGamh(iPs); % plot gam_h
        nfg = nfg+1; figure(nfg);
        pdeplot(points,[],triangles,'XYData',lsf{1},'XYStyle','interp',...
         'ZData',lsf{1},'ZStyle','continuous',...
         'ColorBar','off','Title','Numerical solution'); 
    end
end


%% ========================================================
% EXACT SOLUTION
% =========================================================
if eS{2} % user wanna plot exact solution
    figure(nfg); nfg = nfg+1;
    if eS{3}==2 % 2D
        pdeplot(points,edges,triangles,'XYData',eS{1},'Title','Exact solution');
        plotMesh(eS{4},msh); % plot with mesh
        plotGamh(iPs); % plot gam_h
    elseif eS{3}==3 % 3D
        pdeplot(points,[],triangles,'XYData',eS{1},'XYStyle','interp',...
         'ZData',eS{1},'ZStyle','continuous',...
         'ColorBar','off','Title','Exact solution'); 
        plotMesh(eS{4},msh); % plot with mesh
    else % both 1D and 2D
        pdeplot(points,edges,triangles,'XYData',eS{1},'Title','Exact solution');
        plotMesh(eS{4},msh); % plot with mesh
        plotGamh(iPs); % plot gam_h
        figure(nfg); nfg = nfg+1;
        pdeplot(points,[],triangles,'XYData',eS{1},'XYStyle','interp',...
         'ZData',eS{1},'ZStyle','continuous',...
         'ColorBar','off','Title','Exact solution'); 
    end
end


%% ========================================================
% NUMERICAL SOLUTION
% =========================================================
if nS{2} % user wanna plot numerical solution
    figure(nfg);
    if nS{3}==2 % 2D
        pdeplot(points,edges,triangles,'XYData',nS{1},'Title','Numerical solution');
        plotMesh(nS{4},msh); % plot with mesh
        plotGamh(iPs); % plot gam_h
    elseif nS{3}==3 % 3D
        pdeplot(points,[],triangles,'XYData',nS{1},'XYStyle','interp',...
         'ZData',nS{1},'ZStyle','continuous',...
         'ColorBar','off','Title','Numerical solution'); 
        plotMesh(nS{4},msh); % plot with mesh
    else % both 1D and 2D
        pdeplot(points,edges,triangles,'XYData',nS{1},'Title','Numerical solution');
        plotMesh(nS{4},msh); % plot with mesh
        plotGamh(iPs); % plot gam_h
        nfg = nfg+1; figure(nfg);
        pdeplot(points,[],triangles,'XYData',nS{1},'XYStyle','interp',...
         'ZData',nS{1},'ZStyle','continuous',...
         'ColorBar','off','Title','Numerical solution'); 
    end
end

end