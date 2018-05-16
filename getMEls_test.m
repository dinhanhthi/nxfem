function E = getMEls_test(msh,pa,vold,del)
% TEST: use velocity velo instead of grad of vold
% cf. main_level_simple.m
% get matrix E_ij for level set equation like in Arnold p.221
% E_ij = sum_T of ( xi_j,xi_i + del*dot(gradv,grad xi_i) )_L2
% this technique is diff from finding Global matrix separatelt because we
% can use it (again) in finding load vector
% Input: - vold in Vh std,
%        - del (Arnold Book p.222)
% Output: sparse matrix E


[i1,j1,v1] = getTriplePPTs(msh,pa); % xi*xi
[i2,j2,v2] = getTriplePGGTs_test(vold,msh,pa,del); % xi*del*(gradv*gradxi)

ii = [i1;i2];
jj = [j1;j2];
vv = [v1;v2];

E = sparse(ii,jj,vv);

end