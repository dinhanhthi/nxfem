function plotNXFEMex(nS,eS,mshs,lsf,msh)
% Containing all plots
% This file is used only for export plot to figures (without plotting in a
%               seperated windows)
% Input: - numerical solution & exact solution in standard FEM
%        - 
% Output: figures depending on what users wanna see

nfg = 1; % number of figure
triangles=msh.t; points=msh.p; edges=msh.e;


% ---------------------------------------------------------------
    function plotMesh(yOn,msh) % plot solution with mesh or not?
       if yOn
           hold on
           pdemesh(msh.p,msh.e,msh.t); % plot mesh
           hold off 
       end
    end
% ---------------------------------------------------------------


%% ========================================================
% MESH
% =========================================================
if mshs{1}
    f=figure('visible','off');
%     figure(nfg); nfg = nfg+1;
    pdemesh(points,edges,triangles,'NodeLabels',mshs{2},'ElementLabels',mshs{3});
    print -djpeg results\onlyMesh.jpg;
    close(f);
end


%% ========================================================
% LEVEL SET FUNCTION
% =========================================================
if lsf{2} % user wanna plot level set function
%     figure(nfg); nfg = nfg+1;
    f=figure('visible','off');
    if lsf{3}==2 % 2D
        pdeplot(points,edges,triangles,'XYData',lsf{1},'Title','Level Set function');
        plotMesh(lsf{4},msh); % plot with mesh
        print -djpeg results\levelSetFunction.jpg;
        close(f);
    elseif lsf{3}==3 % 3D
        pdeplot(points,[],triangles,'XYData',lsf{1},'XYStyle','interp',...
         'ZData',lsf{1},'ZStyle','continuous',...
         'ColorBar','off','Title','Level set func');  
        plotMesh(lsf{4},msh); % plot with mesh
        print -djpeg results\levelSetFunction.jpg;
        close(f);
    else % both 1D and 2D
        pdeplot(points,edges,triangles,'XYData',lsf{1},'Title','Numerical solution');
        plotMesh(lsf{4},msh); % plot with mesh
        nfg = nfg+1; figure(nfg);
        pdeplot(points,[],triangles,'XYData',lsf{1},'XYStyle','interp',...
         'ZData',lsf{1},'ZStyle','continuous',...
         'ColorBar','off','Title','Numerical solution'); 
        print -djpeg results\levelSetFunction.jpg;
        close(f);
    end
end


%% ========================================================
% EXACT SOLUTION
% =========================================================
if eS{2} % user wanna plot exact solution
%     figure(nfg); nfg = nfg+1;
    f=figure('visible','off');
    if eS{3}==2 % 2D
        pdeplot(points,edges,triangles,'XYData',eS{1},'Title','Exact solution');
        plotMesh(eS{4},msh); % plot with mesh
        print -djpeg results\exactSolution.jpg;
        close(f);
    elseif eS{3}==3 % 3D
        pdeplot(points,[],triangles,'XYData',eS{1},'XYStyle','interp',...
         'ZData',eS{1},'ZStyle','continuous',...
         'ColorBar','off','Title','Exact solution'); 
        plotMesh(eS{4},msh); % plot with mesh
        print -djpeg results\exactSolution.jpg;
        close(f);
    else % both 1D and 2D
        pdeplot(points,edges,triangles,'XYData',eS{1},'Title','Exact solution');
        plotMesh(eS{4},msh); % plot with mesh
        figure(nfg); nfg = nfg+1;
        pdeplot(points,[],triangles,'XYData',eS{1},'XYStyle','interp',...
         'ZData',eS{1},'ZStyle','continuous',...
         'ColorBar','off','Title','Exact solution'); 
        print -djpeg results\exactSolution.jpg;
        close(f);
    end
end


%% ========================================================
% NUMERICAL SOLUTION
% =========================================================
if nS{2} % user wanna plot numerical solution
%     figure(nfg);
    f=figure('visible','off');
    if nS{3}==2 % 2D
        pdeplot(points,edges,triangles,'XYData',nS{1},'Title','Numerical solution');
        plotMesh(nS{4},msh); % plot with mesh
        print -djpeg results\numericalSolution.jpg;
        close(f);
    elseif nS{3}==3 % 3D
        pdeplot(points,[],triangles,'XYData',nS{1},'XYStyle','interp',...
         'ZData',nS{1},'ZStyle','continuous',...
         'ColorBar','off','Title','Numerical solution'); 
        plotMesh(nS{4},msh); % plot with mesh
        print -djpeg results\numericalSolution.jpg;
        close(f);
    else % both 1D and 2D
        pdeplot(points,edges,triangles,'XYData',nS{1},'Title','Numerical solution');
        plotMesh(nS{4},msh); % plot with mesh
        nfg = nfg+1; figure(nfg);
        pdeplot(points,[],triangles,'XYData',nS{1},'XYStyle','interp',...
         'ZData',nS{1},'ZStyle','continuous',...
         'ColorBar','off','Title','Numerical solution');
        print -djpeg results\numericalSolution.jpg;
        close(f);
    end
end

end