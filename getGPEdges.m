function [eGP,neighborTris] = getGPEdges(tris,phi,msh,pa)
% Find the edges around the interface to compute the penalty terms
% NOTE: we can use this file not only for cut triangles!!!!
% Status: checked with mdl=4, tH=3, regu=0, nSeg=5
% State: check with mdl=3, pa.reguMesh = 0, nSeg=4
% Input: cut triangles + phi + type of tris
% Output: - vector contains ghost-penalty edges (5 x nEs)
%           + lines 1-2: 2 endpoints 
%           + lines 3-4: 2 adjacent triangles (idx in neighborTris)
%           + line 5: type of edge (in which part of Omg?)
%           + lines 6-7: the order of 2 endpoints in the triangle indicated
%               in line 3 (1,2 or 3)
%           + lines 8-9: the order of 2 endpoints in the triangle indicated
%               in line 4 (1,2 or 3)

triangles=msh.t;

%% Find all neighbor triangles of cut triangles
% it includes cut-triangles themselves
%----------------------------------------------------
neighborTris = pdeent(triangles,tris(5,:)); % row array, idx in msh.t
nNBTris = size(neighborTris,2); % number of neighbors of tris


%% find all edges of neighbor triangles
%----------------------------------------------------
eNBtris = zeros(9,3*nNBTris);
eNBtris(1:2,1:nNBTris) = triangles([1,2],neighborTris); % edge i-j
eNBtris(1:2,nNBTris+1:2*nNBTris) = triangles([2,3],neighborTris); % edge j-k
eNBtris(1:2,2*nNBTris+1:3*nNBTris) = triangles([3,1],neighborTris); % edge k-i


%% remember also the 'father' triangle
% idx is in neighborCTs
%----------------------------------------------------
eNBtris(3,1:nNBTris) = 1:nNBTris; % numbering in neighborTris
eNBtris(3,nNBTris+1:2*nNBTris) = 1:nNBTris; % numbering in neighborTris
eNBtris(3,2*nNBTris+1:3*nNBTris) = 1:nNBTris; % numbering in neighborTris



%% rememember position of this vertex in the triangle
%----------------------------------------------------
eNBtris(6:7,1:nNBTris) = repmat([1;2],1,nNBTris); % edge i-j
eNBtris(6:7,nNBTris+1:2*nNBTris) = repmat([2;3],1,nNBTris); % edge j-k
eNBtris(6:7,2*nNBTris+1:3*nNBTris) = repmat([3;1],1,nNBTris); % edge k-i


%% find 2 adjacent triangles
%----------------------------------------------------
tmp = sort(eNBtris(1:2,:),1); % same edges have the same order endpoints
[~,iF,~] = unique(tmp','rows','first','legacy'); % first occurrence
[~,iL,~] = unique(tmp','rows','last','legacy'); % last occurrence
interFL = intersect(iF,iL,'stable'); % not-duplicate positions
posF = setdiff(iF,interFL,'stable');% first occurence of duplicate triangle
posL = setdiff(iL,interFL,'stable'); % last occurence of duplicate triangle

eGP = eNBtris(:,posF); % contain triangle K
eGP(4,:) = eNBtris(3,posL); % contain triangle K', numbering in neighborTris
eGP(8:9,:) = eNBtris([7,6],posL); % position of vertices in K'
% note: i-j in K but j-i in K', that's why we take [7,6] not [6,7]


%% Filter the wrong edges
% eGP must contain all edges belonging to CTs, but there are some edges
% doesn't.
%--------------------------------------------------------
tmp = ismember(neighborTris(eGP(3:4,:))',tris(5,:));
eGP=eGP(:,tmp(:,1)|tmp(:,2));


%% find type of edge
%----------------------------------------------------
neGP = size(eGP,2); % number of ghost-penalty edges
% 1 if phi<0,0 if phi=0, 2 if phi>0
markedPoint = 1*((phi<0)&(abs(phi)>pa.tol))+2*((phi>0)&(abs(phi)>pa.tol)); 
sumEndPointseCTs = sum(markedPoint(eGP(1:2,:))); % row-array
prodEndPointseCTs = prod(markedPoint(eGP(1:2,:))); % row-array

% add line: type of edges of cut triangles
for i=1:neGP
   switch sumEndPointseCTs(i)
       case 1
           eGP(5,i) = 1; % in Omg1
       case 2
           if prodEndPointseCTs(i)==0
               eGP(5,i) = 2; % in Omg2
           else
               eGP(5,i) = 1; % in Omg1
           end
       case 3
           eGP(5,i) = 3; % cut edges
       case 4
           eGP(5,i) = 2; % in Omg2
   end
end

end