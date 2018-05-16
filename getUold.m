function uold = getUold(u,msh)
% Find uold in each subdomain to be used in finding u
% This file firstly used in the file getGMuLinda (get global matrix for u)
% Input: - u in nxfem space
%        - msh to find needed nodes
% Output: - u in standard fem for NCTs1 (don't care at nodes in Omg2)
%         - u in standard fem for NCTs2 (don't care at nodes in Omg1)
%         - u in standard fem for CTs1 (only consider nodes in CTs region)
%         - u in standard fem for CTs2 (only consider nodes in CTs region)
%   nstd x 1

%-------------------------------------------------------------------------
uold.omg1 = sparse(msh.nStd,1); uold.omg2 = sparse(msh.nStd,1);
uold.ct1 = sparse(msh.nStd,1); uold.ct2 = sparse(msh.nStd,1);
% uold.omg1 = zeros(msh.nStd,1); uold.omg2 = zeros(msh.nStd,1);
% uold.ct1 = zeros(msh.nStd,1); uold.ct2 = zeros(msh.nStd,1);

%-------------------------------------------------------------------------
uold.omg1(msh.node.omg1) = u(msh.node.omg1);
%   all nodes in omg1 including nodes on gam
%   Note that, there may be nodes on Gam but not in NCTs1. However, we
%   don't care about them because we only consider the NCTs which don't
%   contain them.

%-------------------------------------------------------------------------
uold.omg2(msh.node.omg2.notCT) = u(msh.node.omg2.notCT);
%   all nodes in omg2 (incl. on gam) but not in cut triangles
uold.omg2(msh.node.CT.iomg2) = u(msh.newNodes(msh.node.CT.iomg2));
%   nodes in CTs and inside omg2
uold.omg2(msh.node.CT.onG) = u(msh.newNodes(msh.node.CT.onG));
%   nodes in CTs and on gam
%   Note that, there may be nodes on Gam but not in NCTs 2. However, we
%   don't care about them because we only consider the NCTs which don't
%   contain them.

%-------------------------------------------------------------------------
uold.ct1(msh.node.CT.all) = u(msh.node.CT.all);
%   all nodes in CTs

%-------------------------------------------------------------------------
uold.ct2(msh.node.CT.all) = u(msh.newNodes(msh.node.CT.all));
%   all nodes in CTs

end