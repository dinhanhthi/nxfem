function uNCT = getUNCT(CTs,typeCT,phi,iPs,pa)
% Find unit normal of all Gam_h of all cut triagles
% State: - checked but there are some round-off errors from matlab
%        - 31/10/18: checked again the algorithm, good, nothing changes
% Input: cut triangles + type of CT + phi + intersection points
% Output: all unit normal vector: matrix of 2 coordinates x nCT

nCT = size(CTs,2); % number of cut triangles
uNCT = zeros(2,nCT);
tol=pa.tol;

    function check = negative(num,tole)
        check = 0;
        if (num<0)&&(abs(num)>tole)
            check=1; 
        end
    end
    function check = positve(num,tole)
        check = 0;
        if (num>0)&&(abs(num)>tole)
            check=1; 
        end
    end

for t=1:nCT
    % type 3, interface passes one vertex
    if typeCT(t)==0
        if negative(phi(CTs(1,t)),tol)||positve(phi(CTs(3,t)),tol)
            % vertex(1) in Omg1 OR vertex(3) in Omg2
           uNCT(:,t) = getUnitNV(iPs(:,2,t),iPs(:,1,t)); % BA
        else
           uNCT(:,t) = getUnitNV(iPs(:,1,t),iPs(:,2,t)); % AB
        end
    % type 1 or 2, interface cuts 2 edges
    elseif negative(phi(CTs(1,t)),tol) % vertex(1) in Omg1
       uNCT(:,t) = getUnitNV(iPs(:,2,t),iPs(:,1,t)); % BA 
    else % vertex(1) in Omg2
       uNCT(:,t) = getUnitNV(iPs(:,1,t),iPs(:,2,t)); % AB
    end
end % end of for nCT

end