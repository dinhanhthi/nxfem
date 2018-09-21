function [ii,jj,vv] = getTriplePPCTs(msh,pa,CTs,CT,sub,varargin)
% Get triplets of int_CTs(K*phi_i*phi_j) where P has form f(w(x,y))*g(u(x,y))*h(x,y)
% Rewrite to  getTriplePPCTs, used for flexible function K
% Old file: getLCTs
% Note: 
%   - it's different from getTriplePPCTs, this function only find 1 vv at once,
%   it's neccessary to indicate sub in the input
%   - max 2 variables w and u
% Status: Checked with the old getTriplePPCTs in getGMgPP and getGMuNewton 
% Related: 
% Input: - cut triangles CTs
%        - CT's info
%        - the choice of form of f,g,h (function handles)
%           + ALL FUNCTION HANDLES MUST BE DEFINED BE4 PARSING TO THIS FUNCTION
%           + h(x,y,pa)
%           + f(w), g(u)
%        - u,w are in stdFEM
%        - sub: 1 (for ij), 2 (for kikj)
% Output: ii (nodes); vv (values) column arrays


%% Prerequisites
points = msh.p;
defaultVar = ones(msh.nStd,1);
iPs=CT.iPs; typeCTs=CT.type;nodeCTs=CT.nodes;
nCTs = size(CTs,2); % number of all not cut triangles in Omg1
ii = zeros(9*nCTs,1); jj = zeros(9*nCTs,1); vv = zeros(9*nCTs,1);


%% Setting up inputParser
p = inputParser; % initial inputParser
addRequired(p,'msh'); % mesh
addRequired(p,'pa'); % fixed parameters
addRequired(p,'CTs'); % not-cut triangles
addRequired(p,'CT'); % not-cut triangles
addRequired(p,'sub'); % domain indication, 1 or 2

dG = @findDefG; % default function g=1
addParameter(p,'gu',dG,@(f) isa(f,'function_handle')); % function handle g of var u (u in NXFEM)
addParameter(p,'u',defaultVar); % u in stdFEM

dF = @findDefF; % default function f=1
addParameter(p,'fw',dF,@(f) isa(f,'function_handle')); % function handle f of var w (w in NXFEM)
addParameter(p,'w',defaultVar); % w in stdFEM

dH = @findDefH; % default function h = 1tri
addParameter(p,'h',dH,@(f) isa(f,'function_handle')); % free function handle

parse(p,msh,pa,CTs,CT,sub,varargin{:});
rP = p.Results;


%% Get triplet
idx = 1;
for t=1:nCTs
    triangle = CTs(:,t); % info of triangle t-th
    iP1 = iPs(:,1,t); % 1st intersection point
    iP2 = iPs(:,2,t); % 2nd intersection point
    typeT = typeCTs(t);
    if typeT==2
        rV = points(:,nodeCTs.eachOmg2(1,t)); 
    else
        rV = points(:,nodeCTs.eachOmg1(1,t)); 
    end
    for i=1:3
        for j=1:3
            ii(idx) = CTs(i,t);
            jj(idx) = CTs(j,t);
            vWhole = getTriplePPWhole(msh,pa,triangle,j,i,...
                'u',rP.u,'gu',rP.gu,'w',rP.w,'fw',rP.fw,'h',rP.h);
            vPart = getTriplePPPart(msh,pa,triangle,j,i,iP1,iP2,rV,...
                    'u',rP.u,'gu',rP.gu,'w',rP.w,'fw',rP.fw,'h',rP.h);
            if ( (sub==1)&&(typeT==2) ) || ( (sub==2)&&(typeT~=2) )
                vv(idx) = vWhole - vPart;
            else
                vv(idx) = vPart;
            end
            idx = idx+1;
        end
    end
end


%% Default function handles
function val = findDefH(xx,yy,pa)
    val= 1;
end

function val = findDefG(uu)
    val= 1;
end

function val = findDefF(ww)
    val= 1;
end


end % the end 
