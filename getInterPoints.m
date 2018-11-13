function interPoints = getInterPoints(cutTriangles,phi,msh,pa)
% Find intersection points between phi and cut triangles
% Note that, interPoints DOES contain point on the interface
% Input: cut triangles, points and phi
% Output: - intersection points : 2 coordinates x 2 cut points x nCT

nCT = size(cutTriangles,2); % number of cut triangles
temp = zeros(2,2,nCT); % (2 coordinates x 2 cut points x nCT) matrix of NaN
cP = zeros(1,nCT); % count the cut points
points=msh.p; tol=pa.tol;

%% we need 2 following functions to avoid the round-off problem in matlab
    function check = equal0(num1,num2,tole)
        if (abs(num1)<tole)||(abs(num2)<tole)
            check=1;
        else
            check=0;
        end
    end
    function check = negative(num1,num2,tole)
        if (num1*num2<0)&&(abs(num1)>tole)&&(abs(num2)>tole)
            check=1; 
        else
            check=0;
        end
    end

%%
for i=1:nCT
    phiTri = phi(cutTriangles(1:3,i)); % phi on 3 vertices (1x3)
    pt(:,:) = points(:,cutTriangles(1:3,i)); % 2 x number of vertices
    if negative(phiTri(1),phiTri(2),tol) % edge 1
        cP(i) = cP(i)+1;
        temp(:,cP(i),i) = getInterPointOnEdge(pt(:,1),pt(:,2),phiTri(1),phiTri(2));
    end
    if negative(phiTri(2),phiTri(3),tol)||equal0(phiTri(2),phiTri(3),tol) % edge 2
        cP(i) = cP(i)+1;
        temp(:,cP(i),i) = getInterPointOnEdge(pt(:,2),pt(:,3),phiTri(2),phiTri(3));
    end
    if negative(phiTri(3),phiTri(1),tol)||equal0(phiTri(3),phiTri(1),tol) % edge 3
        cP(i) = cP(i)+1;
        temp(:,cP(i),i) = getInterPointOnEdge(pt(:,3),pt(:,1),phiTri(3),phiTri(1));
    end
end

interPoints = temp;
end