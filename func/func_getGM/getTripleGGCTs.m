function [ii,jj,v1,v2] = getTripleGGCTs(CTs,areaChildCTs,k1,k2,msh)
% Get triplets of int_CTs(kk*gradPhi_i*gradPhi_j)
% Input: - information on CTs
%        - coefficients kk1,kk2 on each subdomain
% Output: triplet for both cases in each subdomain on CTs (column-arrays)

points=msh.p;

nCTs = size(CTs,2); % number of cut triangles
ii = zeros(9*nCTs,1); % column-arrays
jj = zeros(9*nCTs,1); % column-arrays
v1 = zeros(9*nCTs,1); % column-arrays
v2 = zeros(9*nCTs,1); % column-arrays

% normal vectors
% ve is a 2 coordinates x nCT x 3 edges matrix
ve(:,:,3) = points(:,CTs(2,:)) - points(:,CTs(1,:));
ve(:,:,1) = points(:,CTs(3,:)) - points(:,CTs(2,:));
ve(:,:,2) = points(:,CTs(1,:)) - points(:,CTs(3,:));
areaCTs = 0.5*abs(ve(1,:,3).*ve(2,:,2) - ve(1,:,2).*ve(2,:,3)); 

idx = 0;
for i=1:3
    for j=1:3
        % for basis locating on Omg1 
        ii(idx+1:idx+nCTs) = CTs(i,:);
        jj(idx+1:idx+nCTs) = CTs(j,:);
        v1(idx+1:idx+nCTs) = k1*areaChildCTs(1,:)./(4*areaCTs.^2)...
                .*dot(ve(:,:,j),ve(:,:,i)); % A_ij=A(j,i)
        v2(idx+1:idx+nCTs) = k2*areaChildCTs(2,:)./(4*areaCTs.^2)...
                .*dot(ve(:,:,j),ve(:,:,i)); % A_ij=A(j,i)
        idx = idx+nCTs;
    end % end for j
end % end for i

end