function val = getNormL2oGstd(msh,pa,CTs,CT,err)
% Find ||err_h||_{L^2(Gam_h)} for stdFEM err
% First used in main_level_simple to find norm of signed distance function,
%   cf. Arnold book p.224
% Suppose that CTs is not empty
% Input: - cut triangles' info
%        - err in stdFEM: 1 x nTs
% Output: scalar value

L = []; % take the default value, L = 1 for all cases
[iPP,jPP,vPP1,~,~,~] = getTriplePPoG(CTs,CT.iPs,msh,pa,L);

mA = sparse(iPP,jPP,vPP1,msh.nStd,msh.nStd);

val = err*mA*err';
val = sqrt(val);

end
