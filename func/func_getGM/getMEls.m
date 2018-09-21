function E = getMEls(msh,pa,velo,delT,coef)
% Velocity is optional, del is optional
% Get matrix E_ij for level set equation like in Arnold p.221
% E_ij = coef * sum_T of ( xi_j,xi_i + del*dot(velo,grad xi_i) )_T
% State: checked with the old one in commit f5c821
% This technique is diff from finding Global matrix separately because we
%   can use it (again) in finding load vector
% Related: main_chopp2007.m, main_levelset_simple, getMHls.m, note 6
% Input: - velo(x,y,pa,sub) function handle, sub indicates velo.x or y
%        - delT (Arnold Book p.222) or something else? : 1xnTs
% Output: sparse matrix E to be used in finding matrix of level set eqn



%% Get info
Ts = msh.t; % all triandelNewgles of the mesh
defVx = @(x,y,pa) velo(x,y,1);
defVy = @(x,y,pa) velo(x,y,2);


%% int_Ts phi_j*phi_i
K = []; % take the default K
[i1,j1,v1] = getTriplePPNCTs(msh,pa,Ts,K);


% sum_T of ( xi_j,del*dot(velo,grad xi_i) )_T
delT = delT*coef;
[i2,j2,v2] = getTripleVeloGPP(msh,pa,Ts,delT,'xi',defVx);
[i3,j3,v3] = getTripleVeloGPP(msh,pa,Ts,delT,'yi',defVy);


ii = [i1;i2;i3];
jj = [j1;j2;j3];
vv = [v1;v2;v3];


E = sparse(ii,jj,vv);


end
