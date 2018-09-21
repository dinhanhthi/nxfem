function [ii,ff1,ff2] = getfPhiCTs(msh,pa,CTs,varargin)
% NOT FINISHED YET!!!!
% int_CTs P*phi on CTs where P has form f(w(x,y))*g(u(x,y))*h(x,y)
% Rewrite to  getLCTs, used for flexible function P
% Old file: getLCTs
% Note: 
%   - it's different from getLNCTs where u,w are structure with u.ct1,u.ct2
%   and h depends on sub (to indicate omg1 or omg2)
%   - max 2 variables w and u
% Status: 
% Related: test_getfPhiCTs.m
% Input: - cut triangles CTs
%        - the choice of form of f,g,h (function handles)
%           + ALL FUNCTION HANDLES MUST BE DEFINED BE4 PARSING TO THIS FUNCTION
%           + h(x,y,pa,sub) <-- there is 'sub'
%           + f(w), g(u)
%        - u (u.ct1,u.ct2) and w (w.ct1,w.ct2): all *.ct* is in stdFEM
% Output: ii (nodes); ffi for basis whose support in Omgi (values at nodes), column arrays

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

dG = @findDefG; % default function g=1
addParameter(p,'gu',dG,@(f) isa(f,'function_handle')); % function handle g of var u (u in NXFEM)
addParameter(p,'u',defaultVar); % u in stdFEM

dF = @findDefF; % default function f=1
addParameter(p,'fw',dF,@(f) isa(f,'function_handle')); % function handle f of var w (w in NXFEM)
addParameter(p,'w',defaultVar); % w in stdFEM

dH = @findDefH; % default function h = 1
addParameter(p,'h',dH,@(f) isa(f,'function_handle')); % free function handle
parse(p,msh,pa,CTs,varargin{:});
rP = p.Results;

idx=1;
for t=1:nCTs
    iP1 = iPs(:,1,t); % 1st intersection point
    iP2 = iPs(:,2,t); % 2nd intersection point
    triangle = CTs(:,t);
    for i=1:3 % 3 vertices
        ii(idx) = CTs(i,t);
        if typeCTs(t)==2 % 1 node in Omg2, 2 nodes in Omg1
            rV = points(:,nodeCTs.eachOmg2(1,t)); % the only vertex in Omg2
            Fwhole1 = getLWhole(tri,i,msh,pa,PPW1);
            Fwhole1 = getfPhiWhole(msh,pa,triangle,i,...
                        'u',rP.u,'gu',rP.gu,'w',rP.w,'fw',rP.fw,'h',rP.h);
            Fpart1 = getLPart(tri,i,iP1,iP2,rV,msh,pa,PPP1);
            Fpart1 = getfPhiPart(msh,pa,triangle,i,iP1,iP2,rV,...
                        'u',rP.u,'gu',rP.gu,'w',rP.w,'fw',rP.fw,'h',rP.h);
            ff(idx) = Fwhole1 - Fpart1;
            ff2(idx) = getLPart(tri,i,iP1,iP2,rV,msh,pa,PPP2);
            idx = idx+1;
        else % typeCT = 0 or 4
            rV = points(:,nodeCTs.eachOmg1(1,t)); % the only vertex in Omg1
            Fwhole2 = getLWhole(tri,i,msh,pa,PPW2);
            Fpart2 = getLPart(tri,i,iP1,iP2,rV,msh,pa,PPP2);
            ff2(idx) = Fwhole2 - Fpart2;
            ff(idx) = getLPart(tri,i,iP1,iP2,rV,msh,pa,PPP1);
            idx = idx+1;
        end % end if typeCT
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