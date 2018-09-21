function Fp = getfPhiPart(msh,pa,triangle,vertex,iP1,iP2,rV,varargin)
% int_T P*phi_j on a "part" triangle where P has form f(w(x,y))*g(u(x,y))*h(x,y)
% Rewrite to getLPart, this file can be used for a flexible function P
%   P can be: a function handle, a FE variable or couple of them
% Note: max 2 variables w and u
% Status: check with getfPhiCTs.m
% Related: getfPhiCTs.m, test_getfPhiPart.m. getfPhiWhole.m
% Input: - which node i?
%        - which triangle? (with points to get vertices)
%        - 2 intersection points: iP1, iP2
%        - only 1 remaining vertex: rV
%        - the choice of form of f,g,h (function handles)
%           + ALL FUNCTION HANDLES MUST BE DEFINED BE4 PARSING TO THIS FUNCTION
%           + h(x,y,pa)
%           + f(w), g(u)
%        - u (and w) in stdFEM
% Output: scalar value


points=msh.p;
defaultVar = ones(msh.nStd,1);
% pa.degN: Gaussian quadrature points (for complicated function)
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D
[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points


%% =======================================================================
% Setting up inputParser
%=========================================================================
p = inputParser; % initial inputParser
addRequired(p,'msh'); % mesh
addRequired(p,'pa'); % fixed parameters
addRequired(p,'triangle'); % considered triangle
addRequired(p,'vertex',@(s) ismember(s,[1,2,3])); % considered vertex of triangle t (1, 2 or 3)
addRequired(p,'iP1'); % intersection 1
addRequired(p,'iP2'); % intersection 2
addRequired(p,'rV'); % remaining vertex

dG = @findDefG; % default function g=1
addParameter(p,'gu',dG,@(f) isa(f,'function_handle')); % function handle g of var u (u in NXFEM)
addParameter(p,'u',defaultVar); % u in stdFEM

dF = @findDefF; % default function f=1
addParameter(p,'fw',dF,@(f) isa(f,'function_handle')); % function handle f of var w (w in NXFEM)
addParameter(p,'w',defaultVar); % w in stdFEM

dH = @findDefH; % default function h = 1
addParameter(p,'h',dH,@(f) isa(f,'function_handle')); % free function handle
parse(p,msh,pa,triangle,vertex,iP1,iP2,rV,varargin{:});
rP = p.Results;


%% =======================================================================
% Get fPhiPart
%=========================================================================

v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3
areaT = getAreaTri(v1,v2,v3); % area of triangle

% get u,w in each triangle
uS = rP.u(triangle(1:3));
wS = rP.w(triangle(1:3));

[xiP1h,yiP1h] = getCoorRef(rP.iP1,v1,v2,v3); % cut point 1 in ref coor
[xiP2h,yiP2h] = getCoorRef(rP.iP2,v1,v2,v3); % cut point 2 in ref coor
[xRvh,yRvh] = getCoorRef(rP.rV,v1,v2,v3); % remaining vertex in ref coor
% ref-triangle's part's area
areaTHp = getAreaTri([xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);

Fp = 0;
for k=1:nwt
    % point in intermediate coordinate wrt the Gaussian point in ref coor
    [xHk,yHk] = getCoorSTD(pt(:,k),[xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);
    % shape function N_i at quadrature point in intermediate coordinate
    [shFu,~,~] = getP1shapes(xHk,yHk);
    vuS = uS(1)*shFu(1)+uS(2)*shFu(2)+uS(3)*shFu(3);
    vwS = wS(1)*shFu(1)+wS(2)*shFu(2)+wS(3)*shFu(3);
    [xk,yk] = getCoorSTD([xHk,yHk],v1,v2,v3); % gaussian points in origin part of triangle
    Fp = Fp + rP.fw(vwS)*rP.gu(vuS)*rP.h(xk,yk,pa)*2*areaT*areaTHp*wt(k)*shFu(vertex);
end


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


end