function [idxOnGamCT,idxOmg1CT,idxOmg2CT] = getIdxCT(CTs,phi,pa)
% Get information about the indices of cut triangles
% Input: cut triangles and phi
% Output: - idx of vertices of each triangle in each subdomain
%         - idx of vertices of each triangle on the interface
%     Note that, just store the nodes in EACH triangle, 
%     DON'T care about their positions

nCT = size(CTs,2); % number of cut triangles
phi = phi'; % to be sure tempi is always column vector!
phiCT = phi(CTs(1:3,:)); % phi on vertices of each triangle
cutTriangleVT = CTs(1:3,:); % consider only 3 vertices

% temp0 = find(phiCT==0); % on interface
temp0 = find(abs(phiCT)<pa.tol); % on interface
temp1 = find((phiCT<0)&(abs(phiCT)>pa.tol)); % in Omg1
temp2 = find((phiCT>0)&(abs(phiCT)>pa.tol)); % in Omg2

    function idxOmg = findIdx(tempi,cTvt,nCT,type)
        ni = size(tempi,1);
        idx = zeros(1,nCT)+1;
        tmp = NaN(type,nCT);
        for i=1:ni
            pos = floor((tempi(i)-1)/3)+1; % which colum in CTs (which CT)?
            tmp(idx(pos),pos) = cTvt(tempi(i));
            idx(pos)=idx(pos)+1; % idx(pos)=which CT
        end
        idxOmg=tmp;
    end % end local function

if ~isempty(temp0)
    idx0 = findIdx(temp0,cutTriangleVT,nCT,1); % only 1 point on interface
else
    idx0 = []; % there is no nodes on the interface
end
idx1 = findIdx(temp1,cutTriangleVT,nCT,2); % max 2 points in Omg1
idx2 = findIdx(temp2,cutTriangleVT,nCT,2); % max 2 points in Omg2

% get nodes on each subdomain
idxOnGamCT = idx0; % 1 point x nCT
idxOmg1CT = idx1; % max 2 points x nCT
idxOmg2CT = idx2; % max 2 points x nCT

end % end function