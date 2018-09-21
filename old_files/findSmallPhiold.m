function [tri2del,t2Omg1,t2Omg2] = findSmallPhiold(msh,CTs,iPs,nodeCTs,...
                                        pa,typeCTs,phi)
% Find vertices whose support on the cut triangle is very small
% state: 
% please see (20) in arnold 2008
% Input: - max diam of each triangle
%        - cut triangles CTs, 
%        - nodes infor in CTs: nodeCTs, 
%        - type of CTs, area of CTs
% Output: - vector contains indices of all satisfied triangles IN OLD CTs!!!
%           (these indices are in the OLD CTs array!)
%         - cut triangles transform from CTs to NCTs (t2Omg1, t2Omg2):
%         indices in OLG CTs

points = msh.p;
tH = pa.tH; % \tilde{c}
hT = msh.hT;
nCTs = size(CTs,2); % number of cut triangles
t2Omg1 = zeros(1,nCTs); % initial set of required triangles
t2Omg2 = zeros(1,nCTs);

idx = 0;
for t=1:nCTs % consider all cut triangles
    triangle = CTs(:,t);
    h = hT(CTs(5,t)); % hT of triangle t
    iP1 = iPs(:,1,t); % 1st intersection point
    iP2 = iPs(:,2,t); % 2nd intersection point
    for i=1:3 % consider all 3 vertices of each triangle
        if typeCTs(t)==2 % 2 nodes in Omg1, 1 node in Omg2
            rV = points(:,nodeCTs.eachOmg2(1,t)); % the only vertex in Omg2
            normPhiOmg2 = getNormPhiPart(triangle,i,iP1,iP2,rV,pa,msh);
            normWhole = getNormPhiWhole(triangle,i,pa,msh); % ||.||^2
            normPhiOmg1 = normWhole - normPhiOmg2; % ||.||^2
        else % typeCT = 0 or 4
            rV = points(:,nodeCTs.eachOmg1(1,t)); % the only vertex in Omg1
            normPhiOmg1 = getNormPhiPart(triangle,i,iP1,iP2,rV,pa,msh);
            normWhole = getNormPhiWhole(triangle,i,pa,msh);
            normPhiOmg2 = normWhole - normPhiOmg1;
        end
        if normPhiOmg1/normWhole<tH^2*h^4 && phi(CTs(i,t))<0 &&...
                        abs(phi(CTs(i,t)))>pa.tol
            idx=idx+1;
            t2Omg2(idx)=t;
            break;
        elseif normPhiOmg2/normWhole<tH^2*h^4 && phi(CTs(i,t))>0 &&...
                        abs(phi(CTs(i,t)))>pa.tol
            idx=idx+1;
            t2Omg1(idx)=t;
            break;
        end
    end % end for 1:3
end

t2Omg1 = t2Omg1(t2Omg1~=0);
t2Omg1 = unique(t2Omg1);
t2Omg2 = t2Omg2(t2Omg2~=0);
t2Omg2 = unique(t2Omg2);
tri2del = [t2Omg1,t2Omg2];
end