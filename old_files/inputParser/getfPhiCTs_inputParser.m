function [ii,ff] = getfPhiCTs(msh,pa,CTs,CT,sub,varargin)
% int_CTs P*phi on CTs where P has form f(w(x,y))*g(u(x,y))*h(x,y)
% Rewrite to  getLCTs, used for flexible function P
% Old file: getLCTs
% Note: 
%   - it's different from getLCTs, this function only find 1 ff at once,
%   it's neccessary to indicate sub in the input
%   - max 2 variables w and u
% Status: checked with main_sys_linda for both Ftw and wgu
% Related: test_getfPhiCTs.m
% Input: - cut triangles CTs
%        - CT's info
%        - the choice of form of f,g,h (function handles)
%           + ALL FUNCTION HANDLES MUST BE DEFINED BE4 PARSING TO THIS FUNCTION
%           + h(x,y,pa)
%           + f(w), g(u)
%        - u,w are in stdFEM
% Output: ii (nodes); ff (values) column arrays

nCTs = size(CTs,2); % number of triangleNC
iPs=CT.iPs; typeCTs=CT.type; nodeCTs=CT.nodes;
points=msh.p;
ii = zeros(3*nCTs,1); ff = zeros(3*nCTs,1);
defaultVar = ones(msh.nStd,1);


%% =======================================================================
% Setting up inputParser
%=========================================================================
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

dH = @findDefH; % default function h = 1
addParameter(p,'h',dH,@(f) isa(f,'function_handle')); % free function handle
parse(p,msh,pa,CTs,CT,sub,varargin{:});
rP = p.Results;

idx=1;
for t=1:nCTs
    iP1 = iPs(:,1,t); % 1st intersection point
    iP2 = iPs(:,2,t); % 2nd intersection point
    triangle = CTs(:,t);
    typeT = typeCTs(t);
    if typeT==2
        rV = points(:,nodeCTs.eachOmg2(1,t)); 
    else
        rV = points(:,nodeCTs.eachOmg1(1,t)); 
    end
    for i=1:3
        ii(idx) = CTs(i,t);
        Fwhole = getfPhiWhole(msh,pa,triangle,i,...
                    'u',rP.u,'gu',rP.gu,'w',rP.w,'fw',rP.fw,'h',rP.h);
        Fpart = getfPhiPart(msh,pa,triangle,i,iP1,iP2,rV,...
                    'u',rP.u,'gu',rP.gu,'w',rP.w,'fw',rP.fw,'h',rP.h);
        if ( (sub==1)&&(typeT==2) ) || ( (sub==2)&&(typeT~=2) )
            ff(idx) = Fwhole - Fpart;
        else
            ff(idx) = Fpart;
        end
        idx = idx+1;
    end
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