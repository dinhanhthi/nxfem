function [kap1,kap2] = getKappa(areaChildCTs,pa)
% Defind kappa1, kappa2 like in Becker's paper
% Input: area of each part af cut triangles
% Output: kappa_i : 1 x nCTs

kk1=pa.kk1; kk2=pa.kk2;
nCTs = size(areaChildCTs,2);

if ~pa.kap1per2 % kap1per2=0
%     kap1 = (kk2*areaChildCTs(1,:))...
%                 ./(kk2*areaChildCTs(1,:) + kk1*areaChildCTs(2,:));
% %     kap1 = (areaChildCTs(1,:))...
% %                 ./(areaChildCTs(1,:) + areaChildCTs(2,:)); % hansbo's idea
%     kap2 = 1-kap1;

    %% ghost penalty
    kap1 = zeros(1,nCTs) + kk2/(kk1+kk2);
    kap2 = zeros(1,nCTs) + kk1/(kk1+kk2);
else
    % without Becker's guide
    kap1 = zeros(1,nCTs) + 0.5;
    kap2 = zeros(1,nCTs) + 0.5;
end

end