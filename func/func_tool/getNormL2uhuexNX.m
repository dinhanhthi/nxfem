function val = getNormL2uhuexNX(msh,pa,tris,CT,uh,uex,defPhi)
% Find ||fh-f||_L2(\Omg) in NXFEM
% Related file: getNormL2fhf.m (the same function but in stdFEM)
% First used in main_nxfem (Hansbo's errors between uex and uh)
% Can be used for finding error between discrete and continuous functions
%   in NXFEM
% Input: uh : NXFEM (1 x ndof)
%        uex: function handle (x,y,pa,sub)
%        defPhi(x,y,pa): function handle
% ---------------------------------------------
% Important (3/11): this file is diff (but not diff) from getNormL2fhfSTD, the
%   only diff is on the parser parameters uex(x,y,pa,sub) whereas f(x,y,pa).
%   In function uex, there is uex1, uex2 and we don't know how to determine
%   xk,yk in which subdomain?
% ---------------------------------------------


norm_fh = getNormL2fhNX(uh,tris,CT,msh,pa);
norm_f = getNormL2uexSTD(msh,pa,uex,defPhi); % The only diff is here!!!


ffh = getLf(msh,pa,tris,CT,uex); % ndof x 1
fh_f = dot(uh',ffh);

val = norm_fh^2 + norm_f^2 - 2*fh_f;
val = sqrt(val);

end