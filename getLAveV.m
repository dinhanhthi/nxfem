function L = getLAveV(cpV)
% Get L in ||{gran u}||_{-1/2} = \sum_K hK \int_GamK {gran u}^2
% L in this case: hK*{}{}
% There are relations between kappa and kk. They belong to the same var
% Related file: getNormAveV.m
% IMPORTANT: because we wanna use getTriplePPoG (phi*phi) but there are
% signs minus (-) in lines 3 and 4 whereas in this norm, all signs are +.
% Thus we need to change sign lines 3, 4 of L before applying to this
% function.
% Output: matrix L: 4 x nCTs
%           row 1: Aij
%           row 2: Akikj
%           row 3: Akij
%           row 4: Aikj

kap1 = cpV.kap1; kap2 = cpV.kap2;
nCTs = size(kap1,2);
L = zeros(4,nCTs);

L(1,:) = kap1.^2; % Aij
L(2,:) = kap2.^2; % Akikj
L(3,:) = -kap1.*kap2; % Akij
L(4,:) = -kap1.*kap2; % Aikj

end