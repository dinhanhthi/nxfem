function [ii,jj,vv1,vv2] = getMatEstEtaSCE(CTs,msh,pa,phi)
% Get triplet on cut edges for finding matrix etaS
% Input: given in the input parameters
% Output: triplet to build the final matrix

triangles = msh.t; points = msh.p;

%% Get information about edges
[eCTs,neighborCTs] = getGPEdges(CTs,phi,msh,pa);
icE = find(eCTs(5,:)==3); % consider just the cut edges (this is just the idx)
ne = size(icE,2); % number of considered edges
nbCTs = triangles(:,neighborCTs); % neighbor of cut triangles
gradPhinbCTs = getGradPhi(nbCTs,msh); % grad of phi wrt vertices of nbNCTs
% gradPhinbCTs = 2 coordinates x 3 vertices x nbCTs
kk1 = pa.kk1; kk2 = pa.kk2;


%% setting up
ii = zeros(4*(ne),1); % column-array
jj = zeros(4*(ne),1); % column-array
vv1 = zeros(4*(ne),1); % column-array, for ij
vv2 = zeros(4*(ne),1); % column-array, for k(i)k(j)

%% find triplet
idx=1;
for e=1:ne
    edge = icE(e); % considered edge in eGPs
    vphi1 = phi(eCTs(1,edge)); % phi's value at endpoint 1
    vphi2 = phi(eCTs(2,edge)); % phi's value at endpoint 2
    if (vphi1<0)&&(abs(vphi1)>pa.tol) % enpoint 1 is in Omg1
        eP1 = points(:,eCTs(1,edge)); % endpoint 1
        eP2 = points(:,eCTs(2,edge)); % endpoint 2
    else % endpoint 1 is in Omg2
        eP1 = points(:,eCTs(2,edge)); % endpoint 1
        eP2 = points(:,eCTs(1,edge)); % endpoint 2
        % permute vphi1 and vphi2
        vphi1 = vphi1+vphi2;
        vphi2 = vphi1-vphi2;
        vphi1 = vphi1-vphi2;
    end
    iP = getInterPointOnEdge(eP1,eP2,vphi1,vphi2); % intersection point
    normalVT = getUnitNV(eP1,eP2); % unit normal vector to e
    lenEdge1 = sqrt((eP1(1)-iP(1))^2 +(eP1(2)-iP(2))^2); % length of edge in Omg1
    lenEdge2 = sqrt((eP2(1)-iP(1))^2 +(eP2(2)-iP(2))^2); % length of edge in Omg2
    for i=1:2
        for j=1:2
           ii(idx) = eCTs(i,edge);
           jj(idx) = eCTs(j,edge);
           % [grad_n phi_i]_e
%            jumpGradnPhii = getGradnPhi(eCTs(i+5,edge),neighborCTs(eCTs(3,edge)),eP1,eP2,msh)...
%                             - getGradnPhi(eCTs(i+7,edge),neighborCTs(eCTs(4,edge)),eP1,eP2,msh);
           jumpGradnPhii = ...
               dot(gradPhinbCTs(:,eCTs(i+5,edge),eCTs(3,edge)),normalVT)...
             - dot(gradPhinbCTs(:,eCTs(i+7,edge),eCTs(4,edge)),normalVT); 
            % [grad_n phi_j]_e
%            jumpGradnPhij = getGradnPhi(eCTs(j+5,edge),neighborCTs(eCTs(3,edge)),eP1,eP2,msh)...
%                             - getGradnPhi(eCTs(j+7,edge),neighborCTs(eCTs(4,edge)),eP1,eP2,msh);
           jumpGradnPhij = ...
               dot(gradPhinbCTs(:,eCTs(j+5,edge),eCTs(3,edge)),normalVT)...
             - dot(gradPhinbCTs(:,eCTs(j+7,edge),eCTs(4,edge)),normalVT);  
           vv1(idx) = kk1*lenEdge1*jumpGradnPhij*jumpGradnPhii;
           vv2(idx) = kk2*lenEdge2*jumpGradnPhij*jumpGradnPhii;
           idx = idx+1;
        end
    end
end

end