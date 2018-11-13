function valGij = getGij_WholeTri(areaT,ii,jj,typeG,dim,deg)
% Find integrate of g(Phi_i,Phi_j) on the single whole triangle
% Input: - which phi_i (ii) & phi_j (jj)
%        - what is g(.,.)? (kk)
%        - dimension dim and degree of quadrature
% Output: the value of integrate on this triangle wrt nodes ii and jj 


[wt,pt] = getGaussQuad(dim,deg); % quadrature weights and points
nwt = size(wt,2); % number of Gaussian points

valGij = 0;
for k=1:nwt
    [shapeFunction,~,~] = getP1shapes(pt(1,k),pt(2,k)); % shape functions at quadrature points
    valGij = valGij + areaT*wt(k)*defG(shapeFunction(ii),shapeFunction(jj),typeG);
end

end