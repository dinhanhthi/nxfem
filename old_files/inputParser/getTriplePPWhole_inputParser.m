function val = getTriplePPWhole(msh,pa,triangle,i,j,varargin)
% Find int_T(K*phi_j*phi_i) on whole triangle
%   where K has form f(w(x,y))*g(u(x,y))*h(x,y)
% Rewrite to getTriplePPWhole, this file can be used for a flexible 
%    function K (function handle or FE var or couple of them)
% Note: - max 2 FE variables w and u
%       - rewrite the second times (1st is to use K)
%       - the same idea to getfPhiWhole
% Related: getTriplePPNCTs, getTriplePPCTs
% State: checked along with getTriplePPNCTs2 in test_getTriplePPNCTs
% Input: - which node i,j (1,2 or 3)
%        - which triangle
%        - the choice of form of f,g,h (function handles)
%           + ALL FUNCTION HANDLES MUST BE DEFINED BE4 PARSING TO THIS FUNCTION
%           + h(x,y,pa)
%           + f(w), g(u)
%        - u (and w) in stdFEM
% Output: - the value on whole triangle


%% Prerequisites
points = msh.p;
defaultVar = ones(msh.nStd,1);
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D (polynomial func)
[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points


%% Setting up inputParser
p = inputParser; % initial inputParser
addRequired(p,'msh'); % mesh
addRequired(p,'pa'); % fixed parameters
addRequired(p,'triangle'); % considered triangle
addRequired(p,'i',@(s) ismember(s,[1,2,3])); % phi_i
addRequired(p,'j',@(s) ismember(s,[1,2,3])); % phi_j

dG = @findDefG; % default function g=1
addParameter(p,'gu',dG,@(f) isa(f,'function_handle')); % function handle g of var u (u in NXFEM)
addParameter(p,'u',defaultVar); % u in stdFEM

dF = @findDefF; % default function f=1
addParameter(p,'fw',dF,@(f) isa(f,'function_handle')); % function handle f of var w (w in NXFEM)
addParameter(p,'w',defaultVar); % w in stdFEM

dH = @findDefH; % default function h = 1
addParameter(p,'h',dH,@(f) isa(f,'function_handle')); % free function handle
parse(p,msh,pa,triangle,i,j,varargin{:});
rP = p.Results;


%% Get info
v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3
areaT = getAreaTri(v1,v2,v3); % area of triangle

% get u,w in each triangle
uS = rP.u(triangle(1:3));
wS = rP.w(triangle(1:3));


%% Get int value
val=0;
for k=1:nwt
    [shFu,~,~] = getP1shapes(pt(1,k),pt(2,k)); % N_i at quadrature points
    [xk,yk] = getCoorSTD(pt(:,k),v1,v2,v3);
    vuS = uS(1)*shFu(1)+uS(2)*shFu(2)+uS(3)*shFu(3);
    vwS = wS(1)*shFu(1)+wS(2)*shFu(2)+wS(3)*shFu(3);
    val = val + rP.fw(vwS)*rP.gu(vuS)*rP.h(xk,yk,pa)*areaT*wt(k)*shFu(i)*shFu(j);
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
