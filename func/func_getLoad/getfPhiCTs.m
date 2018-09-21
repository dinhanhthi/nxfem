function [ii,ff1,ff2] = getfPhiCTs(msh,pa,CTs,CT,P)
% int_Omg P*phi on CTs (load vector)
% This file is used for general P (the size of P depends on number of
%       Gaussian points (pa.degN)
% Note: rewriten after 1st rewrite with inputParser (worst performance)
% State: - checked with getLf
%        - (old) alike getLoadCTs, diff at 16th digit after ","
% Input: - cut triangles CTs
%        - CTs' info CT
%        - P.CTW1,P.CTW2 for the whole triangle case
%        - P.CTP1,P.CTP2 for the part triangle case
%        (all are nCTs x nwt)
% Output: ii (nodes); ff1, ff2 (values at nodes). column arrays

iPs=CT.iPs; typeCTs=CT.type; nodeCTs=CT.nodes; points=msh.p;
nCTs = size(CTs,2); % number of cut triangles

% default P
if ~isempty(P)
    PW1=P.CTW1; PW2=P.CTW2; PP1=P.CTP1; PP2=P.CTP2;
else
    dim=2; deg=pa.degN; 
    [wt,~] = getGaussQuad(dim,deg); 
    nwt = size(wt,2);
    Pdf = ones(nCTs,nwt);
    PW1=Pdf.CTW1; PW2=Pdf.CTW2; PP1=Pdf.CTP1; PP2=Pdf.CTP2;
end


ii = zeros(3*nCTs,1); ff1 = zeros(3*nCTs,1); ff2 = zeros(3*nCTs,1);

idx=1;
for t=1:nCTs
    iP1 = iPs(:,1,t); % 1st intersection point
    iP2 = iPs(:,2,t); % 2nd intersection point
    triangle = CTs(:,t);
    PPW1=PW1(t,:); PPW2=PW2(t,:); PPP1=PP1(t,:); PPP2=PP2(t,:);
    for i=1:3 % 3 vertices
        ii(idx) = CTs(i,t);
        if typeCTs(t)==2 % 1 node in Omg2, 2 nodes in Omg1
            rV = points(:,nodeCTs.eachOmg2(1,t)); % the only vertex in Omg2
            Fwhole1 = getfPhiWhole(msh,pa,triangle,i,PPW1);
            Fpart1 = getfPhiPart(msh,pa,triangle,i,iP1,iP2,rV,PPP1);
            ff1(idx) = Fwhole1 - Fpart1;
            ff2(idx) = getfPhiPart(msh,pa,triangle,i,iP1,iP2,rV,PPP2);
            idx = idx+1;
        else % typeCT = 0 or 4
            rV = points(:,nodeCTs.eachOmg1(1,t)); % the only vertex in Omg1
            Fwhole2 = getfPhiWhole(msh,pa,triangle,i,PPW2);
            Fpart2 = getfPhiPart(msh,pa,triangle,i,iP1,iP2,rV,PPP2);
            ff2(idx) = Fwhole2 - Fpart2;
            ff1(idx) = getfPhiPart(msh,pa,triangle,i,iP1,iP2,rV,PPP1);
            idx = idx+1;
        end % end if typeCT
    end % end for vertices
end % end for nCT

end
