function phiphi = intPhiPhi(pointA,pointB,ii,jj,v,pa)
% Find the value of int_pointA^pointB of Phi_iPhi_j
% Input: - intersection points: pointA and pointB (1x2)
%        - ii of Phi_i, jj of Phi_j
%        - vertices v of this triangle: 2 coordinates x 3 vertices
% Output: the value of integral for these i and j.

lenAB = sqrt((pointB(1)-pointA(1))^2+(pointB(2)-pointA(2))^2);

% coordinate of point A in the "reference" triangle
[xHa,yHa] = getCoorRef(pointA,v(:,1),v(:,2),v(:,3));
% coordinate of point B in the "reference" triangle
[xHb,yHb] = getCoorRef(pointB,v(:,1),v(:,2),v(:,3));

% Get Gaussian quadrature weights and points
% pa.degP1D: Gaussian quadrature points in 1D (for polynomial function)
dim=1; deg=pa.degP1D; % quadrature's info
[wt,pt] = getGaussQuad(dim,deg); 
nwt = size(wt,2); % number of Gaussian points

tmp = 0; % temporary variable
for k=1:nwt
    xHaT = xHa + (xHb-xHa)*((1+pt(k))/2); % xHa(t), t=(1+xi)/2
    yHaT = yHa + (yHb-yHa)*((1+pt(k))/2); % yHa(t), t=(1+xi)/2
    [shFu,~,~] = getP1shapes(xHaT,yHaT);
    tmp = tmp + wt(k)*shFu(ii)*shFu(jj); % shape function jj, jj in 1:3
end

phiphi = (1/2)*lenAB*tmp;
end