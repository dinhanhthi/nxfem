function fh = model_test_mshdist
    fh.phiRight = @findPhiRight;
    fh.domain = @findDomain;
    fh.phiWrong = @findPhiWrong;
end

function GeoDom = findDomain()
% Output: a matrix containing all information of domain
    xDomVal = [0 1 1 0]; % x values of points constructing Omega
    yDomVal = [0 0 1 1]; % corresponding y value
    RectDom = [3,4,xDomVal,yDomVal]'; % rectangular domain "3" with "4" sides
    GeoDom = decsg(RectDom);
end

function valPhi=findPhiWrong(xx,yy,pa)
% Define level set function phi
% Input: coordinate of points
% Output: value of phi at points
    pa.xi = 0.15;
    valPhi = (xx-0.5).^2 + yy.^2 - pa.xi^2;
end

function valPhi=findPhiRight(xx,yy,pa)
% Define level set function phi
% Input: coordinate of points
% Output: value of phi at points
    pa.xi = 0.15;
    valPhi = sqrt((xx-0.5).^2 + yy.^2) - pa.xi;
end