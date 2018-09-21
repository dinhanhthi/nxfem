function [ii,jj,vv] = getTriplePPNCTs(msh,pa,NCTs,varargin)
% Get triplets of int_NCTs(K*phi_i*phi_j) where K has form f(w(x,y))*g(u(x,y))*h(x,y)
% Rewrite to getTriplePPNCTs, this file can be used for a flexible 
%    function K (function handle or FE var or couple of them)
% Note: - max 2 FE variables w and u
%       - rewrite the second times (1st is to use K)
%       - the same idea to getfPhiNCTs
% Related: test_getTriplePPNCTs, getGMgPP, getGMuNewton.m
% State: Checked with the old getTriplePPNCTs in getGMgPP and getGMuNewton
% Input: - not-cut triangles NCTs
%        - the choice of form of f,g,h (function handles)
%           + ALL FUNCTION HANDLES MUST BE DEFINED BE4 PARSING TO THIS FUNCTION
%           + h(x,y,pa)
%           + f(w), g(u)
%        - u (and w) in stdFEM
% Output: triplet on NCTs wrt each subdomain (column-arrays)


%% Prerequisites
nNCTs = size(NCTs,2); % number of all not cut triangles in Omg1
ii = zeros(9*nNCTs,1); jj = zeros(9*nNCTs,1); vv = zeros(9*nNCTs,1);
defaultVar = ones(msh.nStd,1);


%% Setting up inputParser
p = inputParser; % initial inputParser
addRequired(p,'msh'); % mesh
addRequired(p,'pa'); % fixed parameters
addRequired(p,'NCTs'); % not-cut triangles

dG = @findDefG; % default function g=1
addParameter(p,'gu',dG,@(f) isa(f,'function_handle')); % function handle g of var u (u in NXFEM)
addParameter(p,'u',defaultVar); % u in stdFEM

dF = @findDefF; % default function f=1
addParameter(p,'fw',dF,@(f) isa(f,'function_handle')); % function handle f of var w (w in NXFEM)
addParameter(p,'w',defaultVar); % w in stdFEM

dH = @findDefH; % default function h = 1
addParameter(p,'h',dH,@(f) isa(f,'function_handle')); % free function handle
parse(p,msh,pa,NCTs,varargin{:});
rP = p.Results;


%% Get triplet
idx=1;
for t=1:nNCTs
    triangle = NCTs(:,t);
    for i=1:3
        for j=1:3
            ii(idx) = NCTs(i,t);
            jj(idx) = NCTs(j,t);
            vv(idx) = getTriplePPWhole(msh,pa,triangle,j,i,...
                    'u',rP.u,'gu',rP.gu,'w',rP.w,'fw',rP.fw,'h',rP.h);
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
