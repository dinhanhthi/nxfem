function H = getMHls_gP(msh,pa,vold,delT,coef)
% NEED TO MODIFY LATER, FOLLOW getGMgPP.m AND getGMuNewton.m
% velocity is grad of Phi
% get matrix H_ij for level set equation like in Arnold p.221
% H_ij = sum_T of ( dot(gradv,grad xi_j),xi_i + del*dot(gradv,grad xi_i) )_L2
% this file computes dt*eps*Hij
% this technique is diff from finding Global matrix separatelt because we
% can use it (again) in finding load vector
% Related: main_chopp2007.m, getMEls.m, note 6
% Input: - vold in Vh std,
%        - del (Arnold Book p.222)
% Output: sparse matrix H

Ts = msh.t; % all triangles of the mesh


%% sum_K of ( del*(grad v*grad xi_j, xi_i) )_K
del1 = ones(1,size(Ts,2)); % del=1
del1 = del1*coef; % include coefficients in level set eqn
[i1,j1,v1] = getTriplePGGTs(vold,msh,pa,del1); 


%% sum_T of ( grad v*grad xi_j,del*(grad v*grad xi_i) )_T
delT = delT*coef; % include coefficients in level set eqn
[i2,j2,v2] = getTripleGGGGTs(vold,msh,delT);

ii = [i1;i2];
jj = [j1;j2];
vv = [v1;v2];

H = sparse(ii,jj,vv);

end
