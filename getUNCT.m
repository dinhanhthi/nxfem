function uNormalsCT = getUNCT(cutTriangles,typeCT,phi,interPoints,pa)
% Find unit normal of all Gam_h of all cut triagles
% State: checked but there are some round-off errors from matlab
% Input: cut triangles + type of CT + phi + intersection points
% Output: all unit normal vector: matrix of 2 coordinates x nCT

nCT = size(cutTriangles,2); % number of cut triangles
tmp = zeros(2,nCT);
tol=pa.tol;

    function check = negative(num,tole)
        if (num<0)&&(abs(num)>tole)
            check=1; 
        else
            check=0;
        end
    end
    function check = positve(num,tole)
        if (num>0)&&(abs(num)>tole)
            check=1; 
        else
            check=0;
        end
    end

for i=1:nCT
    if typeCT(i)==0 % type 3, interface passes one vertex
%         if (phi(cutTriangles(1,i))<0)||(phi(cutTriangles(3,i))>0)
        if negative(phi(cutTriangles(1,i)),tol)||positve(phi(cutTriangles(3,i)),tol)
            % vertex(1) in Omg1 OR vertex(3) in Omg2
           tmp(:,i) = getUnitNV(interPoints(:,2,i),interPoints(:,1,i)); % BA
        else
           tmp(:,i) = getUnitNV(interPoints(:,1,i),interPoints(:,2,i)); % AB
        end
    % type 1 or 2, interface cuts 2 edges
%     elseif phi(cutTriangles(1,i))<0 % vertex(1) in Omg1
    elseif negative(phi(cutTriangles(1,i)),tol)
       tmp(:,i) = getUnitNV(interPoints(:,2,i),interPoints(:,1,i)); % BA 
    else % vertex(1) in Omg2
       tmp(:,i) = getUnitNV(interPoints(:,1,i),interPoints(:,2,i)); % AB
    end
end % end of for nCT

uNormalsCT = tmp;
end