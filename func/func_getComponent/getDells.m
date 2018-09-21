function del = getDells(msh,vold)
% get del coeff in level set equation (Arnold Book p.222)
% del = h*||grad v||_Linf^-1
% don't know exactly, but guess
% norm(grad v,inf) = max_T( max_i( |dx(vi)| + |dy(vi)| ) )
% Realated: main_chopp2007.m, getMEls.m, getMHls.m
% Input: msh and vold in std
% Output: del (scalar)

tris = msh.t;
nTs = size(tris,2);
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

end