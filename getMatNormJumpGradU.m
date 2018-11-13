function mNormJumpGradU = getMatNormJumpGradU(CTs,iPs,msh,pa)
% This function is to compute the matrix of finding value of norm of [grad_n uh] or
% [grad_n u-uh] on Gamma. It's used to check the jump of numerical solution
% Note that the algorithm studies from function getMatrixOnGamCTs
% Input: - cut triangles CTs
%        - intersection points iPs
%        - coefficients coef, it depends on the triangles (eg. getLambda)
% Output: - a matrix


nCTs = size(CTs,2); % number of cut triangles
kk1 = pa.kk1; kk2 = pa.kk2;
newNodes=msh.newNodes;

it = zeros(9*nCTs,1); jt = zeros(9*nCTs,1); % column-arrays
vt1 = zeros(9*nCTs,1); vt2 = zeros(9*nCTs,1); % column-arrays
vt3 = zeros(9*nCTs,1); vt4 = zeros(9*nCTs,1); % column-arrays

idx=1;
for t=1:nCTs
    mt = CTs(5,t); % idx of triangle t in msh.t
    pointA = iPs(:,1,t); % 1st intersection point
    pointB = iPs(:,2,t); % 2nd intersection point
    for i=1:3
        for j=1:3
            it(idx) = CTs(i,t);
            jt(idx) = CTs(j,t);
            % (grad_ni * grad_nj)
            gNgN = getGradnPhi(i,mt,pointA,pointB,msh)*...
                            getGradnPhi(j,mt,pointA,pointB,msh);
            % A_{ij}=a(phi_j,phi_i)
            vt1(idx) = kk1^2*gNgN; % ij
            vt2(idx) = kk2^2*gNgN; % k(i)k(j)
            vt3(idx) = -kk1*kk2*gNgN; % k(i)j
            vt4(idx) = -kk2*kk1*gNgN; % ik(j)
            idx = idx+1;
        end % end for j
    end % end for i
end % end for t

ii = it; jj = jt; vv = vt1;
itmp = newNodes(it); jtmp = newNodes(jt);
ii = [ii;itmp]; jj = [jj;jtmp]; vv = [vv;vt2]; 
it3 = newNodes(it); % column-array
ii = [ii;it3]; jj = [jj;jt]; vv = [vv;vt3];
jt4 = newNodes(jt); % column-array
ii = [ii;it]; jj = [jj;jt4]; vv = [vv;vt4];

mNormJumpGradU = sparse(ii,jj,vv);
end