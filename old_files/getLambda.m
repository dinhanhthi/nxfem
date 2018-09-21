function lambda = getLambda(areaChildCTs,CTs,iPs,pa,msh)
% Defind lambda like in Burman's paper
% Input: - area of each part af cut triangles
%        - cut triangles (CTs)
%        - intersection points: 2 coordinates x 2 cut points x nCT
%        - msh.hT: diam of all triangles 1 x nTs
%        - hTmax: max of diams of all triangles
% Output: lambda: 1 x nCTs

%% setting up
nCTs = size(iPs,3);
hT = msh.hT; hTmax = msh.hTmax;
lenAB = zeros(1,nCTs); % nCTs x 1
lenAB(1,:) = ((iPs(1,2,:)-iPs(1,1,:)).^2+(iPs(2,2,:)-iPs(2,1,:)).^2).^(0.5);
kk1=pa.kk1; kk2=pa.kk2;
lambdaH = pa.lambdaH;
areaCTs = getAreaTris(CTs,msh);

%% just constant
% lambda = zeros(1,size(iPs,3)) + lambdaH; 


%% Barrau's thesis p.31
% it works!
% hTCTs = hT(CTs(5,:)); % consider only on cut triangles
% lambda = (lambdaH*max(kk1,kk2))./(hTCTs);


%% Becker's article
% lambda = lambdaH*max(kk1*areaChildCTs(1,:)./areaCTs,kk2*areaChildCTs(2,:)./areaCTs).*lenAB./areaCTs;
 

%% like in Barrau's thesis (p.32)
% lambda = (lambdaH*kk1*kk2.*lenAB)...
%             ./(kk2*areaChildCTs(1,:) + kk1*areaChildCTs(2,:));


%% Annavarapu's article (mentioned in Burman2014)
% hTCTs = hT(CTs(5,:)); % consider only on cut triangles
% lambda = (lambdaH*(hTCTs/hTmax)*kk1*kk2.*lenAB)...
%             ./(kk2*areaChildCTs(1,:) + kk1*areaChildCTs(2,:));


%% ghost penalty (Burman 2014)
lambda = zeros(1,nCTs) + 4*lambdaH*(kk1*kk2)/(kk1+kk2);
% lambda = zeros(1,nCTs) + 2*(kk1*kk2)/(kk1+kk2);

end