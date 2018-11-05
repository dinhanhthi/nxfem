function E = getMEls_gP(msh,pa,vold,delT,coef)
% Velocity is grad of Phi
% Get matrix E_ij for level set equation like in Arnold p.221
% E_ij = sum_T of ( xi_j,xi_i + del*dot(gradv,grad xi_i) )_T
% This technique is diff from finding Global matrix separatelt because we
%   can use it (again) in finding load vector
% Related: main_chopp2007.m, getMHls.m, note 6
% Input: - vold in Vh std,
%        - delT (Arnold Book p.222) or something else? : 1xnTs
% Output: sparse matrix E
% ----------------------------------------------------------
% Update 26/10/18: This file used to compute grad v on vertices of each
%   triangle. We now use pdegrad to find grad v on the center of each
%   triangle (getMElsgP.m), i.e. grad v is constant on whole triangle and discont at edge.
% ----------------------------------------------------------

Ts = msh.t;

%% int_Ts phi_j*phi_i
K = []; % take the default K
[i1,j1,v1] = getTriplePPNCTs(msh,pa,Ts,K);

%% sum_K of ( xi_j,del*(grad v*grad xi_i) )_K
delT = delT*coef;
[i2,j2,v2] = getTriplePGGTs(vold,msh,pa,delT);

ii = [i1;i2];
jj = [j1;j2];
vv = [v1;v2];

E = sparse(ii,jj,vv);

end
