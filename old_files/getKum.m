function K = getKum(tris,pa)
% get all necessary K for global matrix
% case: int_Omg mu*phi*phi (mu is constant from pa.mu)
% This is K for: mu (constant)
% related file: getGMmPP.m
% Input:
% Ouput: - K.NC1, K.NC2: for NCTs
%        - K.CTW1, K.CTW2, K.CTP1, K.CTP2: for CTs

NCTs1=tris.NCTs1; NCTs2=tris.NCTs2; CTs=tris.CTs;

dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
[wt,~] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2);

K.NC1 = zeros(size(NCTs1,2),nwt) + pa.mu1; % matrix K with all constants mu
K.NC2 = zeros(size(NCTs2,2),nwt) + pa.mu2;

K.CTW1 = zeros(size(CTs,2),nwt) + pa.mu1;
K.CTW2 = zeros(size(CTs,2),nwt) + pa.mu2;
K.CTP1 = K.CTW1;
K.CTP2 = K.CTW2;

end