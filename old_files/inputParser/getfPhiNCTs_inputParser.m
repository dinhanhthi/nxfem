function [ii,ff] = getfPhiNCTs(msh,pa,NCTs,varargin)
% int_NCTs P*phi on NCTs where P has form f(w(x,y))*g(u(x,y))*h(x,y)
% Rewrite to  getLNCTs, used for flexible function P
% Old file: getLNCTs
% Note: max 2 variables w and u
% Status: checked with main_sys_linda Ftw for boths NCTs1 and NCTs2
%               - and for getLfwg including u,w
% Related: test_getfPhiNCTs.m
% Input: - not cut triangles NCTs
%        - the choice of form of f,g,h (function handles)
%           + ALL FUNCTION HANDLES MUST BE DEFINED BE4 PARSING TO THIS FUNCTION
%           + h(x,y,pa)
%           + f(w), g(u)
%        - u (and w) in stdFEM
% Output: ii (nodes), ff (values at nodes), column arrays

nNCTs = size(NCTs,2); % number of triangleNC
ii = zeros(3*nNCTs,1); ff = zeros(3*nNCTs,1);
defaultVar = ones(msh.nStd,1);


%% =======================================================================
% Setting up inputParser
%=========================================================================
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

idx=1;
for t=1:nNCTs
    triangle = NCTs(:,t);
    for i=1:3 % 3 vertices
        ii(idx) = NCTs(i,t);
        ff(idx) = getfPhiWhole(msh,pa,triangle,i,...
                        'u',rP.u,'gu',rP.gu,'w',rP.w,'fw',rP.fw,'h',rP.h);
        idx = idx+1;
    end % end for vertices
end % end for nNCTs


%% =======================================================================
% Default function handle values
%=========================================================================
function val = findDefH(xx,yy,pa)
    val= 1;
end

function val = findDefG(uu)
    val= 1;
end

function val = findDefF(ww)
    val= 1;
end


end % end of main function