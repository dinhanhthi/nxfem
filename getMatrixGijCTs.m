function [ii,jj,vv1,vv2] = getMatrixGijCTs(CTs,iPs,typeCTs,...
                                nodeCTs,areaCTs,typeG,msh,pa)
% get triplet matrix of term int_Omg(G(phi_i,phi_j)) on not cut triangles
% Input: - information of CTs 
%        - type of function G (typeG)
% Output: triplet in CTs wrt each subdomain (column-arrays)

points=msh.p;
dim=2; deg=pa.degP2D; %Gaussian quadrature points in 2D
nCTs = size(CTs,2); % number of all cut triangles
% *2 'cause we compute for both cases in Omg1 and Omg2
ii = zeros(9*nCTs,1); % column-array
jj = zeros(9*nCTs,1); % column-array
vv1 = zeros(9*nCTs,1); % column-array
vv2 = zeros(9*nCTs,1); % column-array

idx=1;
for t=1:nCTs
    iP1 = iPs(:,1,t); % 1st intersection point
    iP2 = iPs(:,2,t); % 2nd intersection point
    triangle = CTs(:,t);
    for i=1:3 % 3 vertices
        for j=1:3 % 3 vertices
            ii(idx) = CTs(i,t);
            jj(idx) = CTs(j,t);
            if typeCTs(t)==2 % 1 node in Omg2, 2 nodes in Omg1
                vInOmg2 = points(:,nodeCTs.eachOmg2(1,t)); % the only node in Omg2
                vv2(idx) = getGij_PartTri(triangle,areaCTs(t),j,i,...
                                typeG,iP1,iP2,vInOmg2,dim,deg,msh); % A_ij=A(i,j)
                vWhole = getGij_WholeTri(areaCTs(t),j,i,typeG,dim,deg); % A_ij=A(i,j)
                vv1(idx) = vWhole - vv2(idx);
                idx = idx+1;
            else % typeCT = 0 or 4
                vInOmg1 = points(:,nodeCTs.eachOmg1(1,t)); % the only vertex in Omg1
                vv1(idx) = getGij_PartTri(triangle,areaCTs(t),j,i,...
                                typeG,iP1,iP2,vInOmg1,dim,deg,msh); % A_ij=A(i,j)
                vWhole = getGij_WholeTri(areaCTs(t),j,i,typeG,dim,deg); % A_ij=A(i,j)
                vv2(idx) = vWhole - vv1(idx);
                idx = idx+1;
            end % end if typeCT
        end % end for j
    end % end for i
end % end for t

end