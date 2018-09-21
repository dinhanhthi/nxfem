function H = getMHls_gP(msh,pa,vold,del,dt,ep)
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
%        - dt time step
%        - ep: theta in theory of full scheme (take ep=1/2)
% Output: sparse matrix H


[i1,j1,v1] = getTriplePGGTs_gP(vold,msh,pa,dt*ep); % dt*(ep*gradv*gradphi)*xi
[i2,j2,v2] = getTripleGGGGTs_gP(vold,msh,del*dt*ep); % dt*ep*del*(gradv*gradphi)*(gradv*gradxi)

ii = [i1;i2];
jj = [j1;j2];
vv = [v1;v2];

H = sparse(ii,jj,vv);

end
