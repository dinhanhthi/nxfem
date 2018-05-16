function CT = getInfoCTs(CTs,phi,msh,pa)
% Find all information about cut triangles
% Input: cut triangles and level set function phi
% Output: Please see details in each section
%           - nodeCTs.each* : nodes on each subdomain/on the interface of 
%                               each triangle in CTs
%           - nodeCTs.onG/Omg1/Omg2 : nodes on each subdomain/on the  
%                               interface of all cut triangles
%             they are all are row arrays
%           - type : 1 x nCT
%           - areaChild: 2 Omg x nCT
%           - iPs: 2 coordinates x 2 cut points x nCT
%           - uN: 2 coordinates x nCT


%% ========================================================
% Get idx of triangles on interface, in Omg1 or in Omg2
% ========================================================= 
% nodes on each subdomain/on gam of EACH triangle in CTs
[idxOnGamCTs,idxOmg1CTs,idxOmg2CTs] = getIdxCT(CTs,phi,pa);
nodeCTs.eachOnG = idxOnGamCTs;  % max 1 point x nCTs, on the interface
nodeCTs.eachOmg1 = idxOmg1CTs; % max 2 points x nCTs, in Omg1
nodeCTs.eachOmg2 = idxOmg2CTs; % max 2 points x nCTs, in Omg2

% nodes on the interface for CTs, row-array
nodeCTs.onG = unique(idxOnGamCTs(~isnan(idxOnGamCTs)));
% nodes in Omg1 for CTs, column-array
nodeCTs.Omg1 = unique(idxOmg1CTs(~isnan(idxOmg1CTs)));
nodeCTs.Omg1 = nodeCTs.Omg1'; % convert to row-array
% nodes in Omg2 for CTs, column-array
nodeCTs.Omg2 = unique(idxOmg2CTs(~isnan(idxOmg2CTs)));
nodeCTs.Omg2 = nodeCTs.Omg2'; % convert to row-array
% number of nodes around the interface for CTs
nodeCTs.n = numel(nodeCTs.onG)+numel(nodeCTs.Omg1)+numel(nodeCTs.Omg2); 
% all nodes in CTs, row-array
nodeCTs.all = [nodeCTs.onG nodeCTs.Omg1 nodeCTs.Omg2];
% there may be some nodes on interface but not in nodeCTs.all (because they
%           are not in the CTs)

CT.nodes = nodeCTs;

%% ========================================================
% Get type of cut triangles
% ========================================================= 
% "0": 1 point on interface; "2": 2 points in Omg1, "4": 2 points in Omg2
% 1 if phi<0,0 if phi=0, 2 if phi>0
% markedPoint = 1*(phi<0)+2*(phi>0);
markedPoint = 1*((phi<0)&(abs(phi)>pa.tol))+2*((phi>0)&(abs(phi)>pa.tol)); 
CT.type = prod(markedPoint(CTs(1:3,:))); % row-array


%% ========================================================
% Get intersection points
% ========================================================= 
CT.iPs = getInterPoints(CTs,phi,msh,pa); 


%% ========================================================
% Get area of each part of cut triangles
% ========================================================= 
CT.areaChild = getAreaChild(CTs,CT.iPs,idxOmg1CTs,idxOmg2CTs,CT.type,msh);


%% ========================================================
% Find unit normal of all Gam_h of all cut triagles
% ========================================================= 
CT.uN = getUNCT(CTs,CT.type,phi,CT.iPs,pa);

end