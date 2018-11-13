function Fti = getLoadWholeTri(triangle,fType,areaT,ii,dim,deg,msh,pa,model)
% Find load vector on the whole triangle
% Input: - which node ii?
%        - which triangle? (with points to get vertices)
%        - triangle's area
%        - which f? (f1 or f2) fType
%        - dimension dim and degree of quadrature
% Output: load vector's value wrt this triangle at node ii

defF = model.defF;
points=msh.p;

[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points
v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3

Fti=0;
for k=1:nwt
    [xk,yk] = getCoorSTD(pt(:,k),v1,v2,v3);
    [shapeFunction,~,~] = getP1shapes(pt(1,k),pt(2,k)); % N_i at quadrature points
    Fti = Fti + areaT*wt(k)*defF(xk,yk,fType,pa)*shapeFunction(ii);
end

end