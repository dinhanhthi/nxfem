function val = getNormL2fhf(msh,pa,fh,f)
% Find ||fh-f||_L2(\Omg)
% First used in finding errors of level set function
% Can be used for finding error between discrete and continuous functions in stdFEM
% Input: fh : stdFEM (1 x nPs)
%        f: function handle (x,y,pa)


norm_fh = getNormL2std(msh,pa,fh);
norm_f = getNormL2fstd(msh,pa,f);


Ts = msh.t; points = msh.p;
nTs = size(Ts,2);
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
[wt,pt] = getGaussQuad(dim,deg);
nwt = size(wt,2);
P = zeros(nTs,nwt);
for t = 1:nTs
    triangle = Ts(:,t);
    v1 = points(:,triangle(1)); % vertex 1
    v2 = points(:,triangle(2)); % vertex 2
    v3 = points(:,triangle(3)); % vertex 3
    for k=1:nwt
        [xk,yk] = getCoorSTD(pt(:,k),v1,v2,v3);
        P(t,k) = f(xk,yk,pa);
    end
end
[ii,ff] = getfPhiNCTs(msh,pa,Ts,P);
ffh = accumarray(ii,ff);
fh_f = dot(fh',ffh);

val = norm_fh^2 + norm_f^2 - 2*fh_f;

val = sqrt(val);

end