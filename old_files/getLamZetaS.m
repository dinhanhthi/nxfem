function lambda = getLamZetaS(areaChildCTs,CTs,iPs,pa,msh)
% Defind coefficients like in Barrau's thesis, p. 32, for the term zeta_S
% Input: - area of each part af cut triangles
%        - cut triangles (CTs)
%        - intersection points: 2 coordinates x 2 cut points x nCT
%        - hT: diam of all triangles 1 x nTs
%        - hTmax: max of diams of all triangles
% Output: lambda: 1 x nCTs

nCTs = size(iPs,3);
lenAB = zeros(1,nCTs); % nCTs x 1
lenAB(1,:) = ((iPs(1,2,:)-iPs(1,1,:)).^2+(iPs(2,2,:)-iPs(2,1,:)).^2).^(0.5);
kk1=pa.kk1; kk2=pa.kk2;
aCTs = getAreaTris(CTs,msh); % 1xnCTs
        
% like in Barrau's thesis (p.32)
lambda = (pa.lambdaH*kk1*kk2.*aCTs)...
            ./(lenAB.*(kk2*areaChildCTs(1,:) + kk1*areaChildCTs(2,:)));

end