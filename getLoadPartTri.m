function [Fk1,Fk2] = getLoadPartTri(triangle,areaT,ii,iP1,iP2,...
                        remainVertex,dim,deg,msh,pa,model)
% Find load vector on a part of triangle (triangle-shape one)
% Input: - which node ii?
%        - which triangle? (with points to get vertices)
%        - triangle's area: areaT
%        - ref-triangle's part's area: areaTp
%        - which f? (f1 or f2): fType
%        - 2 intersection points: iP1, iP2
%        - only 1 remaining vertex: vInOmg1
%        - dimension dim and degree of quadrature
% Output: - load vector wrt f_1: Fk1
%         - load vector wrt f_2: Fk2

defF = model.defF;
points=msh.p;

[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points

v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3

[xiP1h,yiP1h] = getCoorRef(iP1,v1,v2,v3); % cut point 1 in reference coordinate
[xiP2h,yiP2h] = getCoorRef(iP2,v1,v2,v3); % cut point 2 in reference coordinate
[xRvh,yRvh] = getCoorRef(remainVertex,v1,v2,v3); % remaining vertex in ref coor

% ref-triangle's part's area
areaTHp = getAreaTri([xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);

Fk1 = 0; Fk2 = 0;
for k=1:nwt
    % point in intermediate coordinate wrt the Gaussian point in reference coordiante
    [xHk,yHk] = getCoorSTD(pt(:,k),[xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);
    % point in source coordinate wrt the Gaussian point in reference coordiante
    [xk,yk] = getCoorSTD([xHk,yHk],v1,v2,v3);
    % shape function N_i at quadrature point in intermediate coordinate
    [shapeFunction,~,~] = getP1shapes(xHk,yHk);
    Fk1 = Fk1 + 2*areaT*areaTHp*wt(k)*defF(xk,yk,1,pa)*shapeFunction(ii);
    Fk2 = Fk2 + 2*areaT*areaTHp*wt(k)*defF(xk,yk,2,pa)*shapeFunction(ii);
end

end