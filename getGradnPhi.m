function gnp = getGradnPhi(i,t,pA,pB,msh)
% This function is to find grad_n phi of triangle t at vertex i on segment
%           AB
% The order of A and B is very important, unit normal alway points the left
% hand side
% Status: checked!
% Input: - triangle t in msh.t
%        - vertex i
%        - 2 endpoints pA, pB (x,y)
% Output: - value (scalar)

gP = getGrad(i,t,msh); % grad phi
uN = getUnitNV(pA,pB); % unit normal vector

gnp = dot(gP,uN); % grad * n
    
end