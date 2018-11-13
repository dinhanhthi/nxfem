function normWhole = getNormPhiWhole(areaT,ii,dim,deg)
% Find norm of phi_ii on the whole triangle
% state: checked with mdl=4, tH=3, regu=0, nSeg=5
% Input: - which node ii?
%        - triangle's area
%        - dimension dim and degree of quadrature
% Output: - norm of phi_ii on the whole triangle

[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points

normWhole=0;
for k=1:nwt
    [shapeFunction,~,~] = getP1shapes(pt(1,k),pt(2,k)); % N_i at quadrature points
    normWhole = normWhole + areaT*wt(k)*shapeFunction(ii)*shapeFunction(ii);
end

end