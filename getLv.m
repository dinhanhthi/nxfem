function L = getLv(cpV)
% Get L in int_Gam L*phi*phi
% L in this case: lam*{}{} + lam*kap1*kap2*[][]
% Related file: getGMvAA.m, getTriplePPoG.m
% Input:
% Output: matrix L: 4 x nCTs
%           row 1: Aij
%           row 2: Akikj
%           row 3: Akij
%           row 4: Aikj

kap1 = cpV.kap1; kap2 = cpV.kap2;
lambda = cpV.lambda;
nCTs = size(lambda,2);
L = zeros(4,nCTs);

L(1,:) = lambda.*kap1;
L(2,:) = lambda.*kap2;
% L(3,:) = 0 from zeros
% L(4,:) = 0 from zeros

end