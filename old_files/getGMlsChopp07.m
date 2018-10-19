function A = getGMlsChopp07(msh,pa,vold,dt,ep)
% get global matrix for level set function phi in chopp07
% related files: main_chopp2007.m, note 6
% NOTE: this matrix is considered in standard finite element
% Input: - vold in Vh standard
%        - dt time step
%        - ep: theta in theory of full scheme
% Output: global matrix (Vh standard)

tris = msh.t;

% ----------------------------------------
% delta (Arnold Book p.222)
% ----------------------------------------
% don't know exactly, but guess
% norm(grad v,inf) = max_T( max_i( |dx(vi)| + |dy(vi)| ) )
maxGv = 0; % initial max
gP = getGradPhi(tris,msh); % 2 coordinates x 3 vertices x nTris
for t=1:nTs
    v1t = abs(vold(tris(1,t)))*( abs(gP(1,1,t)) + abs(gP(2,1,t)) ); % vertex 1
    v2t = abs(vold(tris(2,t)))*( abs(gP(1,2,t)) + abs(gP(2,2,t)) ); % vertex 2
    v3t = abs(vold(tris(3,t)))*( abs(gP(1,3,t)) + abs(gP(2,3,t)) ); % vertex 3
    % max_i on triangle t
    maxGv = max(maxGv,max(v1t,v2t)); % compare with maxGv on previous t
    maxGv = max(maxGv,v3t); % compare with v3t
end
del = msh.hTmax/maxGv; % del=h*norm(grad v,inf)^{-1}


% ----------------------------------------
% get triples
% ----------------------------------------
[i1,j1,v1] = getTriplePPTs(msh,pa); % phi*xi
[i2,j2,v2] = getTriplePGGTs(vold,msh,pa,del); % phi*del*(gradv*gradxi)
[i3,j3,v3] = getTriplePGGTs(vold,msh,pa,dt*ep); % dt*(ep*gradv*gradphi)*xi
[i4,j4,v4] = getTripleGGGGTs(vold,msh,del*dt*ep); % del*(gradv*gradphi)*(gradv*gradxi)

ii = [i1;i2;i3;i4];
jj = [j1;j2;j3;j4];
vv = [v1;v2;v3;v4];


% ----------------------------------------
% build matrix
% ----------------------------------------
A = sparse(ii,jj,vv);

end