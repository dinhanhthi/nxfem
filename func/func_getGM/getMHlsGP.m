function H = getMHlsGP(msh,pa,gP,delT,coef)
% Velocity is grad of Phi: gradPhi obtained from pdegrad
% This file is diff from getMHls_gP because gP in this case is constant on
%   each triangle (and discont at triangle's edge). The file getMHls_gP
%   considers grapPhi at the vertices.
% Get matrix H_ij for level set equation like in Arnold p.221
% H_ij = sum_T of ( dot(gradv,grad xi_j),xi_i + del*dot(gradv,grad xi_i) )_L2
% This file computes coef*Hij
% This technique is diff from finding Global matrix separatelt because we
%   can use it (again) in finding load vector
% Related: main_chopp2007.m, getMEls.m, note 6
% Input: - gP: gP.x (1 x nTs), gP.y (1 x nTs)
%        - del (Arnold Book p.222)
%        - coef goes with Hij
% Output: sparse matrix H
% ----------------------------------------------------------
% Update 26/10/18: The OLD file used to compute grad v on vertices of each
%   triangle. We now use pdegrad to find grad v on the center of each
%   triangle (this file), i.e. grad v is constant on whole triangle and discont at edge.
% ----------------------------------------------------------

Ts = msh.t; % all triangles of the mesh


%% sum_K of ( del*(grad v*grad xi_j, xi_i) )_K
del1 = ones(1,size(Ts,2)); % del=1
del1 = del1*coef; % include coefficients in level set eqn
[i1,j1,v1] = getTriplePGGTs(gP,msh,pa,del1); 


%% sum_T of ( grad v*grad xi_j,del*(grad v*grad xi_i) )_T
delT = delT*coef; % include coefficients in level set eqn
[i2,j2,v2] = getTripleGGGGTs(gP,msh,delT);

ii = [i1;i2];
jj = [j1;j2];
vv = [v1;v2];

H = sparse(ii,jj,vv);

end