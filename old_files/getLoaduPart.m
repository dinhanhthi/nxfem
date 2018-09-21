function Fp = getLoaduPart(triangle,i,iP1,iP2,rV,dim,deg,kk,...
                                        msh,pa,defF,uoldT,wST)
% Find load vector on a part of triangle (triangle-shape one, for eqn u)
% Input: - which node ii?
%        - which triangle? (with points to get vertices)
%        - ref-triangle's part's area: areaTp
%        - 2 intersection points: iP1, iP2
%        - only 1 remaining vertex: rV
%        - dimension dim and degree of quadrature
%        - uoldT and wST
% Output: load vector's value on a part of triangle at node i

points=msh.p; typeG = pa.typeG;

[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points

v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3
areaT = getAreaTri(v1,v2,v3); % area of triangle

[xiP1h,yiP1h] = getCoorRef(iP1,v1,v2,v3); % cut point 1 in reference coordinate
[xiP2h,yiP2h] = getCoorRef(iP2,v1,v2,v3); % cut point 2 in reference coordinate
[xRvh,yRvh] = getCoorRef(rV,v1,v2,v3); % remaining vertex in ref coor

% ref-triangle's part's area
areaTHp = getAreaTri([xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);

Fp = 0;
for k=1:nwt
    % point in intermediate coordinate wrt the Gaussian point in reference coordiante
    [xHk,yHk] = getCoorSTD(pt(:,k),[xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);
    % point in source coordinate wrt the Gaussian point in reference coordiante
    [xk,yk] = getCoorSTD([xHk,yHk],v1,v2,v3);
    % shape function N_i at quadrature point in intermediate coordinate
    [shFu,~,~] = getP1shapes(xHk,yHk);
    vuold = uoldT(1)*shFu(1)+uoldT(2)*shFu(2)+uoldT(3)*shFu(3);
    guold = defG(vuold,typeG); % g(uold)
    vwS = wST(1)*shFu(1)+wST(2)*shFu(2)+wST(3)*shFu(3); % w
    Fp = Fp + 2*areaT*areaTHp*wt(k)*defF(xk,yk,kk,pa)*shFu(i)...
            - 2*areaT*areaTHp*wt(k)*vwS*guold*shFu(i);
end

end