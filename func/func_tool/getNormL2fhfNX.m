function val = getNormL2fhfNX(msh,pa,tris,CT,fh,f)
% Find ||fh-f||_L2(\Omg) in NXFEM
% Related file: getNormL2fhf.m (the same function but in stdFEM)
% First used in main_nxfem (Hansbo's errors between uex and uh)
% Can be used for finding error between discrete and continuous functions
%   in NXFEM
% Input: fh : NXFEM (1 x ndof)
%        f: function handle (x,y,pa,sub)
% ---------------------------------------------
% Important (3/11): this f in this case is NOT for the exact solution
%   uex(x,y,pa,sub), for this type of function (handle), cf.
%   getNormL2uhuexNX.m
% ---------------------------------------------


norm_fh = getNormL2fhNX(fh,tris,CT,msh,pa);
norm_f = getNormL2fSTD(msh,pa,f);  % The only diff is here!!!


ffh = getLf(msh,pa,tris,CT,f); % ndof x 1
fh_f = dot(fh',ffh);

val = norm_fh^2 + norm_f^2 - 2*fh_f;
val = sqrt(val);

end