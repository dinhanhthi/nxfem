function [ii,jj,vv] = getGhostPenalty(CTs,phi,msh,pa,cp)
% Get 2 ghost penalty terms j_1, j_2
% Status: (29-10-18) big modification on "hT(idxNBTris(eGP(3,edge)))", the
%   old one is "hT(neighborCTs(eGP(3,edge)))" (it's wrong)
% Input: cut triangles + phi + hT of all triangles
% Output: triplet to contribute to build global matrix

hT = msh.hT;

%% Get information about edges
[eGP,idxNBTris] = getGPEdges(CTs,phi,msh,pa);

% classify into 2 sub domains
% idx in eGP
j1 = find(eGP(5,:)==1); % not-cut edges for term j_1
jc = find(eGP(5,:)==3); % cut edges for both terms j_1 & j_2
j2 = find(eGP(5,:)==2); % not-cut edges for term j_2

nej1 = size(j1,2); % number of not-cut edges for j_1
nej2 = size(j2,2); % number of not-cut edges for j_2
nejc = size(jc,2); % number of cut edges
neighborCTs = msh.t(:,idxNBTris); % neighbor cut triangles (4 x number of neighbor triangles)
gradPhinbCTs = getGradPhi(neighborCTs,msh); % grad of phi wrt vertices of neighborCTs
% gradPhinbCTs = 2 coordinates x 3 vertices x number of neighborCTs


%% setting up
ii = zeros(4*(nej1+nej2+2*nejc),1); % column-array
jj = zeros(4*(nej1+nej2+2*nejc),1); % column-array
vv = zeros(4*(nej1+nej2+2*nejc),1); % column-array

kk1 = cp.kk1; kk2 = cp.kk2;
gam1 = pa.gam1; gam2 = pa.gam2; tol=pa.tol;
points = msh.p;
newNodes = msh.newNodes;

%% on cut edges (for both j_1, j_2)
idx=1;
for e=1:nejc
    edge = jc(e); % considered edge in eGPs, edge = idx in eGP
    if (phi(eGP(1,edge))<0)&&(abs(phi(eGP(1,edge)))>tol)
        eP1 = points(:,eGP(1,edge)); % endpoint 1 in Omg1: 2x1
        eP2 = points(:,eGP(2,edge)); % endpoint 2 in Omg2: 2x1
        phi1 = phi(eGP(1,edge)); % phi at ep1
        phi2 = phi(eGP(2,edge)); % phi at ep2
    else
        eP2 = points(:,eGP(1,edge)); % endpoint 2 in Omg2: 2x1
        eP1 = points(:,eGP(2,edge)); % endpoint 1 in Omg1: 2x1
        phi2 = phi(eGP(1,edge)); % phi at ep2
        phi1 = phi(eGP(2,edge)); % phi at ep1
    end
    iP = getiPsonE(eP1,eP2,phi1,phi2); % intersection point
    normalVT = getUnitNV(eP1,eP2); % unit normal vector to e
    len1 = sqrt((eP1(1)-iP(1))^2 +(eP1(2)-iP(2))^2); % length of part of edge in Omg1 (for j1)
    len2 = sqrt((eP2(1)-iP(1))^2 +(eP2(2)-iP(2))^2); % length of part of edge in Omg1 (for j1)
    hE = max( hT(idxNBTris(eGP(3,edge))),hT(idxNBTris(eGP(4,edge))) ); % max h of 2 h in 2 adjacent triangles
    for i=1:2
       for j=1:2
           % [grad_n phi_i]_e
           jumpGradnPhii = ...
               dot(gradPhinbCTs(:,eGP(i+5,edge),eGP(3,edge)),normalVT)...
             - dot(gradPhinbCTs(:,eGP(i+7,edge),eGP(4,edge)),normalVT); 
           % [grad_n phi_j]_e
           jumpGradnPhij = ...
               dot(gradPhinbCTs(:,eGP(j+5,edge),eGP(3,edge)),normalVT)...
             - dot(gradPhinbCTs(:,eGP(j+7,edge),eGP(4,edge)),normalVT);
         
           % for j_1
           ii(idx) = eGP(i,edge); % for j1
           jj(idx) = eGP(j,edge); % for j1
           vv(idx) = kk1*gam1*hE*len1*jumpGradnPhij*jumpGradnPhii;
           idx=idx+1;
           
           % for j_2
           ii(idx) = newNodes(eGP(i,edge)); % for j2
           jj(idx) = newNodes(eGP(j,edge)); % for j2
           vv(idx) = kk2*gam2*hE*len2*jumpGradnPhij*jumpGradnPhii;
           idx=idx+1;
       end
    end
end % nejc

%% on not-cut edges in Omg1 (for j_1)
for e=1:nej1
    edge = j1(e); % considered edge in eGPs, edge = idx in eGP
    eP1 = points(:,eGP(1,edge)); % endpoint 1 in Omg1
    eP2 = points(:,eGP(2,edge)); % endpoint 2 in Omg1
    normalVT = getUnitNV(eP1,eP2); % unit normal vector to e
    lenEdge = sqrt((eP1(1)-eP2(1))^2 +(eP1(2)-eP2(2))^2); % length of edge
    hE = max( hT(idxNBTris(eGP(3,edge))),hT(idxNBTris(eGP(4,edge))) ); % max h of 2 h in 2 adjacent triangles
    for i=1:2
        for j=1:2
           ii(idx) = eGP(i,edge);
           jj(idx) = eGP(j,edge);
           % [grad_n phi_i]_e
           jumpGradnPhii = ...
               dot(gradPhinbCTs(:,eGP(i+5,edge),eGP(3,edge)),normalVT)...
             - dot(gradPhinbCTs(:,eGP(i+7,edge),eGP(4,edge)),normalVT); 
           % [grad_n phi_j]_e
           jumpGradnPhij = ...
               dot(gradPhinbCTs(:,eGP(j+5,edge),eGP(3,edge)),normalVT)...
             - dot(gradPhinbCTs(:,eGP(j+7,edge),eGP(4,edge)),normalVT);  
           vv(idx) = kk1*gam1*hE*lenEdge*jumpGradnPhij*jumpGradnPhii;
           idx = idx+1;
        end
    end
end

%% on not-cut edges in Omg2 (for j_2)
for e=1:nej2
    edge = j2(e); % considered edge in eGPs
    eP1 = points(:,eGP(1,edge)); % endpoint 1 in Omg2
    eP2 = points(:,eGP(2,edge)); % endpoint 2 in Omg2
    normalVT = getUnitNV(eP1,eP2); % unit normal vector to e
    lenEdge = sqrt((eP1(1)-eP2(1))^2 +(eP1(2)-eP2(2))^2); % length of edge
    hE = max( hT(idxNBTris(eGP(3,edge))),hT(idxNBTris(eGP(4,edge))) ); % max h of 2 h in 2 adjacent triangles
    for i=1:2
        for j=1:2
           ii(idx) = newNodes(eGP(i,edge));
           jj(idx) = newNodes(eGP(j,edge));
           % [grad_n phi_i]_e
           jumpGradnPhii = ...
               dot(gradPhinbCTs(:,eGP(i+5,edge),eGP(3,edge)),normalVT)...
             - dot(gradPhinbCTs(:,eGP(i+7,edge),eGP(4,edge)),normalVT); 
           % [grad_n phi_j]_e
           jumpGradnPhij = ...
               dot(gradPhinbCTs(:,eGP(j+5,edge),eGP(3,edge)),normalVT)...
             - dot(gradPhinbCTs(:,eGP(j+7,edge),eGP(4,edge)),normalVT);  
           vv(idx) = kk2*gam2*hE*lenEdge*jumpGradnPhij*jumpGradnPhii;
           idx = idx+1;
        end
    end
end

end