function E = getMElsGP(msh,pa,gP,delT,coef)
% Velocity is grad of Phi: gradPhi obtained from pdegrad
% This file is diff from getMEls_gP because gP in this case is constant on
%   each triangle (and discont at triangle's edge). The file getMEls_gP
%   considers grapPhi at the vertices.
% Get matrix E_ij for level set equation like in Arnold p.221
% E_ij = sum_T of ( xi_j,xi_i + del*dot(gradv,grad xi_i) )_T
% This file actually computes coef*Eij
% This technique is diff from finding Global matrix separately because we
%   can use it (again) in finding load vector
% Related: main_chopp2007.m, getMHls.m, note 6
% Input: - gP: gP.x (1 x nTs), gP.y (1 x nTs) 
%        - delT (Arnold Book p.222) or something else? : 1xnTs
%        - coef goest with Eij
% Output: sparse matrix E
% ----------------------------------------------------------
% Update 26/10/18: The OLD file used to compute grad v on vertices of each
%   triangle. We now use pdegrad to find grad v on the center of each
%   triangle (this file), i.e. grad v is constant on whole triangle and discont at edge.
% ----------------------------------------------------------

Ts = msh.t;

%% int_Ts phi_j*phi_i
K = []; % take the default K
[i1,j1,v1] = getTriplePPNCTs(msh,pa,Ts,K);

%% sum_K of ( xi_j,del*(grad v*grad xi_i) )_K
delT = delT*coef;
[i2,j2,v2] = getTriplePGGTs(gP,msh,pa,delT);

ii = [i1;i2];
jj = [j1;j2];
vv = [v1;v2];

E = sparse(ii,jj,vv);

end
