function Fti = getfPhiWhole(msh,pa,triangle,vertex,varargin)
% int_T P*phi_j on a whole triangle where P has form f(w(x,y))*g(u(x,y))*h(x,y)
% Rewrite to getLWhole, this file can be used for a flexible function P
%   P can be: a function handle, a FE variable or couple of them
% Note: max 2 variables w and u
% Status: - Checked for h, diff from getLWhole, there is no 'sub' in defF
%         - Checked for f,g
% Related: getfPhiNCTs.m, test_getfPhiWhole.m
% Input: - which basis function vertex?
%        - which triangle? (with points to get vertices)
%        - the choice of form of f,g,h (function handles)
%           + ALL FUNCTION HANDLES MUST BE DEFINED BE4 PARSING TO THIS FUNCTION
%           + h(x,y,pa)
%           + f(w), g(u)
%        - u (and w) in stdFEM
% Output: scalar value


points = msh.p; % points of the mesh
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

dG = @findDefG; % default function g=1
addParameter(p,'gu',dG,@(f) isa(f,'function_handle')); % function handle g of var u (u in NXFEM)
addParameter(p,'u',defaultVar); % u in stdFEM

dF = @findDefF; % default function f=1
addParameter(p,'fw',dF,@(f) isa(f,'function_handle')); % function handle f of var w (w in NXFEM)
addParameter(p,'w',defaultVar); % w in stdFEM

dH = @findDefH; % default function h = 1
addParameter(p,'h',dH,@(f) isa(f,'function_handle')); % free function handle
parse(p,msh,pa,triangle,vertex,varargin{:});
rP = p.Results;


%% =======================================================================
% Get fPhiWhole
%=========================================================================

v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3
areaT = getAreaTri(v1,v2,v3); % area of triangle

% get u,w in each triangle
uS = rP.u(triangle(1:3));
wS = rP.w(triangle(1:3));



%% =======================================================================
% Get fPhiWhole
%=========================================================================
Fti=0;
for k=1:nwt
    [shFu,~,~] = getP1shapes(pt(1,k),pt(2,k)); % N_i at quadrature points
    [xk,yk] = getCoorSTD(pt(:,k),v1,v2,v3);
    vuS = uS(1)*shFu(1)+uS(2)*shFu(2)+uS(3)*shFu(3);
    vwS = wS(1)*shFu(1)+wS(2)*shFu(2)+wS(3)*shFu(3);
    Fti = Fti + rP.fw(vwS)*rP.gu(vuS)*rP.h(xk,yk,pa)*areaT*wt(k)*shFu(vertex);
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


end % end of main function