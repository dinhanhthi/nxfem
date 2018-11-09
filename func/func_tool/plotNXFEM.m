function nf = plotNXFEM(msh,pa,phi,iPs,cnf,varargin)
% Plot the mesh or solutions of problems using NXFEM method
% Arguments can be called manually and flexibly
% Status: there's no input checker
% Input:
%   - msh: required arg containing [p,t,e] of the mesh
%   - pa: global parameters
%   - phi: value of phi at nodes
%   - iPs: intersection points (to plot Gam_h)
%   - cnf: current number of plot (if there are alreay other plots)
%   - sol: solutions, it's optional in case you just wanna plot the mesh
%   - `withMesh`: plot sol with the mesh or not: false (default) or true
%   - `withGamh`: plot Gam_h or not? (interface): true (default) or true
%   - `dim`: plot in 2 (2D, default) or 3 (3D)?
%   - `eLabel`: display element label: `off` (default) or `on`
%   - `nLabel`: display node label: `off` (default) or `on`
%   - `title`: title of the plot: manual
%   - `export`: export to a .jpg file with the name the same with title
% Output: - the plot
%         - the current number of figure: nf (in order to do other plots)

%% =======================================================================
% Setting up inputParser
%=========================================================================
points = msh.p; edges = msh.e; triangles = msh.t;
p = inputParser; % initial inputParser
addRequired(p,'msh',@isstruct); % mesh
addRequired(p,'pa',@isstruct); % mesh
addRequired(p,'phi'); % intersection points (to plot Gam_h)
addRequired(p,'iPs'); % intersection points (to plot Gam_h)
addRequired(p,'cnf',@isnumeric); % number of plot
addOptional(p,'sol',[]); % solution to be plotted
addParameter(p,'withMesh',0,@islogical); % plot with mesh or not?
addParameter(p,'withGamh',1,@islogical); % plot with interface or not?
addParameter(p,'dim',2,@(s) ismember(s,[2,3])); % 2D or 3D
addParameter(p,'eleLabel','off',@(s) ismember(s,{'on','off'})); 
    % display element label?
addParameter(p,'nodeLabel','off',@(s) ismember(s,{'on','off'})); 
    % display element label?
addParameter(p,'title','no title'); % plot's title
addParameter(p,'export',false,@islogical); 
    % export to a .jpg file with the name the same with title
parse(p,msh,pa,phi,iPs,cnf,varargin{:});
nf = cnf; % number of plot
resu = p.Results;

%% =======================================================================
% Check and do the plot
%=========================================================================

if isempty(resu.sol) % plot the mesh only
    nf = nf+1; figure(nf);
    resu.withMesh = 1;
else % there is a sol and we will plot it
    nf = nf+1; figure(nf);
    switch resu.dim
        case 2 % 2D
            pdeplot(points,edges,triangles,'XYData',resu.sol,'Title',resu.title);
            if (resu.export) && (~strcmp(resu.title,'no title'))
                print(strcat('results\',resu.title),'-dpng'); % export to file
            end
        case 3 % 3D
            pdeplot(points,[],triangles,'XYData',resu.sol,'XYStyle',...
                'interp','ZData',resu.sol,'ZStyle','continuous',...
                'ColorBar','off','Title',resu.title);
            if (resu.export) && (~strcmp(resu.title,'no title'))
                print(strcat('results\',resu.title),'-dpng'); % export to file
            end
    end
end

% plot the mesh or not
%-------------------------------------------------------------------------
if resu.withMesh
    if isempty(resu.sol)
        title('The mesh');
    end
    hold on
    pdemesh(points,edges,triangles,'NodeLabels',resu.nodeLabel,...
                'ElementLabels',resu.eleLabel);
    if resu.withGamh % plot with Gam_h or not
        
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
        
        hold off
    end
    if (resu.export) && (~strcmp(resu.title,'no title'))
        print(strcat('results\',resu.title),'-dpng'); % export to file
    end
    hold off
end

hold off



end