function L = getLAveGn(cpV,hTCTs)
% Get L in ||{kgran u}||_{-1/2} = \sum_K hK \int_GamK {kgran u}^2
% L in this case: hK*{}{}
% There are relations between kappa and kk. They belong to the same var
% Related file: getNormAveGn.m
% Input: - hTCTs : 1 x nCTs
% Output: matrix L: 4 x nCTs
%           row 1: Aij
%           row 2: Akikj
%           row 3: Akij
%           row 4: Aikj

kap1 = cpV.kap1; kap2 = cpV.kap2;
kk1 = cpV.kk1; kk2 = cpV.kk2;
nCTs = size(hTCTs,2);
L = zeros(4,nCTs);

L(1,:) = hTCTs(1,:).*(kap1.^2).*(kk1.^2); % Aij
L(2,:) = hTCTs(1,:).*(kap2.^2).*(kk2.^2); % Akikj
L(3,:) = hTCTs(1,:).*(kap1.*kap2).*(kk1.*kk2); % Akij
L(4,:) = hTCTs(1,:).*(kap1.*kap2).*(kk1.*kk2); % Aikj

end