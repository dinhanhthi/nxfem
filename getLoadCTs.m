function [FCT1,FCT2] = getLoadCTs(CTs,iPs,areaCTs,typeCTs,nodeCTs,msh,pa,model)
% Get small-local load vector for cut triangles
% State: checked by hand: code follows idea + check with double integral
% Input: - interPoints: 2 coor x 2 cut points x nCTs
%        - CTs:  5 x nCT (3 first line contain vertices)
%        - points: 2 x nP
%        - nodeCTs : nodes on each subdomain/on the interface of each triangle in CTs 
% Output: load vector defined on cut triangles

% pa.degN: Gaussian quadrature points (for complicated function)
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
points=msh.p;

nCTs = size(CTs,2); % number of cut triangles
ii = zeros(3*nCTs,1); ff1 = zeros(3*nCTs,1); ff2 = zeros(3*nCTs,1);

idx=1;
for t=1:nCTs
    iP1 = iPs(:,1,t); % 1st intersection point
    iP2 = iPs(:,2,t); % 2nd intersection point
    triangle = CTs(:,t);
    for i=1:3 % 3 vertices
        ii(idx) = CTs(i,t);
        if typeCTs(t)==2 % 1 node in Omg2, 2 nodes in Omg1
            vInOmg2 = points(:,nodeCTs.eachOmg2(1,t)); % the only vertex in Omg2
            % load vector's value wrt basis whose support in Omg2 (triangle-shape)
            [Fpart,ff2(idx)] = getLoadPartTri(triangle,areaCTs(t),...
                                    i,iP1,iP2,vInOmg2,dim,deg,msh,pa,model);
            % load vector's value wrt basis whose support in Omg1 (quadrilateral-shape)
            Fwhole = getLoadWholeTri(triangle,1,areaCTs(t),i,dim,deg,msh,pa,model);
            ff1(idx) = Fwhole - Fpart;
            idx = idx+1;
        else % typeCT = 0 or 4
            vInOmg1 = points(:,nodeCTs.eachOmg1(1,t)); % the only vertex in Omg1
            % load vector's value wrt basis whose support in Omg1 (triangle-shape)
            [ff1(idx),Fpart] = getLoadPartTri(triangle,areaCTs(t),...
                                    i,iP1,iP2,vInOmg1,dim,deg,msh,pa,model);
            % load vector's value wrt basis whose support in Omg2 (quadrilateral-shape)
            Fwhole = getLoadWholeTri(triangle,2,areaCTs(t),i,dim,deg,msh,pa,model);
            ff2(idx) = Fwhole - Fpart;
            idx = idx+1;
        end % end if typeCT
    end % end for vertices
end % end for nCT

FCT1 = accumarray(ii,ff1);
FCT2 = accumarray(ii,ff2);

end