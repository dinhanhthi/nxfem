function [ii,jj,vv] = getMatEstEtaSNCE(NCTs,sub,phi,msh,pa)
% Get triplet on not-cut edges in each subdomain for finding matrix etaS
% This file is used for not-cut edges in both subdomains
% Input: given in the input parameters
% Output: triplet to build the final matrix

%% Get information about edges
[eNCTs,neighborNCTs] = getGPEdges(NCTs,phi,msh,pa);
% eNCTs = all not-cut edges in Omg_i 
if sub==1 % NCTs in Omg1
    kk = pa.kk1;
else % NCTs in Omg2
    kk = pa.kk2;
end
ne = size(eNCTs,2); % number of considered edges
nbNCTs = msh.t(:,neighborNCTs); % neighbor of not cut triangles
gradPhinbNCTs = getGradPhi(nbNCTs,msh); % grad of phi wrt vertices of nbNCTs
% gradPhinbNCTs = 2 coordinates x 3 vertices x nbNCTs
points = msh.p;

%% setting up
ii = zeros(4*(ne),1); % column-array
jj = zeros(4*(ne),1); % column-array
vv = zeros(4*(ne),1); % column-array


%% find triplet
idx=1;
for edge=1:ne
    eP1 = points(:,eNCTs(1,edge)); % endpoint 1
    eP2 = points(:,eNCTs(2,edge)); % endpoint 2
    normalVT = getUnitNV(eP1,eP2); % unit normal vector to e
    lenEdge = sqrt((eP1(1)-eP2(1))^2 +(eP1(2)-eP2(2))^2); % length of edge
    for i=1:2
        for j=1:2
           ii(idx) = eNCTs(i,edge);
           jj(idx) = eNCTs(j,edge);
           % [grad_n phi_i]_e
           jumpGradnPhii = ...
               dot(gradPhinbNCTs(:,eNCTs(i+5,edge),eNCTs(3,edge)),normalVT)...
             - dot(gradPhinbNCTs(:,eNCTs(i+7,edge),eNCTs(4,edge)),normalVT); 
           % [grad_n phi_j]_e
           jumpGradnPhij = ...
               dot(gradPhinbNCTs(:,eNCTs(j+5,edge),eNCTs(3,edge)),normalVT)...
             - dot(gradPhinbNCTs(:,eNCTs(j+7,edge),eNCTs(4,edge)),normalVT);  
           vv(idx) = kk*lenEdge^2*jumpGradnPhij*jumpGradnPhii;
           idx = idx+1;
        end
    end
end
end