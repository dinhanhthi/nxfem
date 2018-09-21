function H = getMHls(msh,pa,velo,delT,coef)
% Velocity is optional, del is optional
% Get matrix H_ij for level set equation like in Arnold p.221
% H_ij = coef * sum_T of ( dot(velo,grad xi_j),xi_i + del*dot(velo,grad xi_i))_T
% This file computes dt*eps*Hij, not just Hij
% Status: - checked with the old one in commit f5c821
%         - checked with getMHls_const.m
% This technique is diff from finding Global matrix separatelt because we
%       can use it (again) in finding load vector
% Related: main_chopp2007.m, main_levelset_simple, getMEls.m, note 6
% Input: - velo(x,y,pa,sub) function handle, sub indicates velo.x or y
%        - delT (Arnold Book p.222) or something else? : 1xnTs
%        - dt: time step
%        - ep: epsilon in theory of full scheme (take ep=1/2)
% Output: sparse matrix H to be used in finding matrix of level set eqn
%           note that, it contains also dt*epsilon


%% Get info
Ts = msh.t; % all triangles of the mesh
defVx = @(x,y,pa) velo(x,y,1); % velo_x
defVy = @(x,y,pa) velo(x,y,2); % velo_y


%% sum_T of ( dot(velo,grad xi_j),xi_i )_T
del1 = ones(1,size(Ts,2)); % del=1
del1 = del1*coef; % include coefficients in level set eqn
[i1,j1,v1] = getTripleVeloGPP(msh,pa,Ts,del1,'xj',defVx);
[i2,j2,v2] = getTripleVeloGPP(msh,pa,Ts,del1,'yj',defVy);


%% sum_T of ( dot(velo,grad xi_j),del*dot(velo,grad xi_i) )_T
delT = delT*coef; % include coefficients in level set eqn

defVx2 = @(x,y,pa) velo(x,y,1)*velo(x,y,1); % velo_x^2
[i3,j3,v3] = getTripleGGVelo(msh,pa,Ts,delT,'xi','xj',defVx2);

defVy2 = @(x,y,pa) velo(x,y,2)*velo(x,y,2); % velo_y^2
[i4,j4,v4] = getTripleGGVelo(msh,pa,Ts,delT,'yi','yj',defVy2);

defVxVy = @(x,y,pa) velo(x,y,1)*velo(x,y,2); % velo_x*velo_y
[i5,j5,v5] = getTripleGGVelo(msh,pa,Ts,delT,'yi','xj',defVxVy);
[i6,j6,v6] = getTripleGGVelo(msh,pa,Ts,delT,'xi','yj',defVxVy);


%% get triple
ii = [i1;i2;i3;i4;i5;i6];
jj = [j1;j2;j3;j4;j5;j6];
vv = [v1;v2;v3;v4;v5;v6];

H = sparse(ii,jj,vv);

end
