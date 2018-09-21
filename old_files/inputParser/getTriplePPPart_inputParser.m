function val = getTriplePPPart(msh,pa,triangle,i,j,iP1,iP2,rV,varargin)
% Find int_T(K*phi*phi) on part triangle (triangle shape one)
%   where K has form f(w(x,y))*g(u(x,y))*h(x,y)
% Rewrite to getTriplePPPart, this file can be used for a flexible 
%    function K (function handle or FE var or couple of them)
% Note: - max 2 FE variables w and u
%       - rewrite the second times (1st is to use K)
%       - the same idea to getfPhiPart
% Related: getTriplePPCTs
% State: 
% Input: - which triangle (with points to get vertices)
%        - which node i,j
%        - 2 intersection points: iP1, iP2
%        - remaining vertex: rV
%        - the choice of form of f,g,h (function handles)
%           + ALL FUNCTION HANDLES MUST BE DEFINED BE4 PARSING TO THIS FUNCTION
%           + h(x,y,pa)
%           + f(w), g(u)
%        - u (and w) in stdFEM
% Output: the value of integrate on each part of triangle wrt nodes i,j


%% Prerequisites
points = msh.p;
dim=2; deg=pa.degN; % Gaussian quadrature points in 2D (polynomial func)
[wt,pt] = getGaussQuad(dim,deg); % 2D and degree 2 (3 points)
nwt = size(wt,2); % number of Gaussian points
defaultVar = ones(msh.nStd,1);


%% Setting up inputParser
p = inputParser; % initial inputParser
addRequired(p,'msh'); % mesh
addRequired(p,'pa'); % fixed parameters
addRequired(p,'triangle'); % considered triangle
addRequired(p,'i',@(s) ismember(s,[1,2,3])); % phi_i
addRequired(p,'j',@(s) ismember(s,[1,2,3])); % phi_j
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

parse(p,msh,pa,triangle,i,j,iP1,iP2,rV,varargin{:});
rP = p.Results;


%% Get info
v1 = points(:,triangle(1)); % vertex 1
v2 = points(:,triangle(2)); % vertex 2
v3 = points(:,triangle(3)); % vertex 3
areaT = getAreaTri(v1,v2,v3); % area of triangle  

[xiP1h,yiP1h] = getCoorRef(iP1,v1,v2,v3); % iP1 in reference coordinate
[xiP2h,yiP2h] = getCoorRef(iP2,v1,v2,v3); % iP2 in reference coordinate
[xRvh,yRvh] = getCoorRef(rV,v1,v2,v3); % remaining vertex in ref coor

% ref-triangle's part's area
areaTHp = getAreaTri([xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);

% get u,w in each triangle
uS = rP.u(triangle(1:3));
wS = rP.w(triangle(1:3));


%% Get int value
val = 0;
for k=1:nwt
    % point in intermediate coordinate wrt the Gaussian point in ref coor
    [xHk,yHk] = getCoorSTD(pt(:,k),[xiP1h,yiP1h],[xiP2h,yiP2h],[xRvh,yRvh]);
    % shape functions at quadrature point in intermediate coordinate
    [shFu,~,~] = getP1shapes(xHk,yHk);
    vuS = uS(1)*shFu(1)+uS(2)*shFu(2)+uS(3)*shFu(3);
    vwS = wS(1)*shFu(1)+wS(2)*shFu(2)+wS(3)*shFu(3);
    [xk,yk] = getCoorSTD([xHk,yHk],v1,v2,v3); % gaussian points in origin part of triangle
    val = val + rP.fw(vwS)*rP.gu(vuS)*rP.h(xk,yk,pa)*2*areaT*areaTHp*wt(k)*shFu(i)*shFu(j);
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
