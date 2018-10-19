function wS = getWsep(tw,msh,k1,k2)
% Find w in each subdomain to be used in finding u
% Not that w=beta*tw (it's kk*tw)
% Status: 
% We learn the idea from getUold.m
% This file firstly used in the file getGMuLinda (get global matrix for u)
% Input: - u in nxfem space
%        - msh to find needed nodes
%        - k1,k2 are beta
% Output: - u in standard fem for NCTs1 (don't care at nodes in Omg2)
%         - u in standard fem for NCTs2 (don't care at nodes in Omg1)
%         - u in standard fem for CTs1 (only consider nodes in CTs region)
%         - u in standard fem for CTs2 (only consider nodes in CTs region)

%-------------------------------------------------------------------------
wS.omg1 = sparse(msh.nStd,1); % column
wS.omg1(msh.node.omg1) = k1*tw(msh.node.omg1);
%   all nodes in omg1 including nodes on gam
%   Note that, there may be nodes on Gam but not in NCTs1. However, we
%   don't care about them because we only consider the NCTs which don't
%   contain them.

%-------------------------------------------------------------------------
wS.omg2 = sparse(msh.nStd,1);
wS.omg2(msh.node.omg2.notCT) = k2*tw(msh.node.omg2.notCT);
%   all nodes in omg2 (incl. on gam) but not in cut triangles
wS.omg2(msh.node.CT.iomg2) = k2*tw(msh.newNodes(msh.node.CT.iomg2));
%   nodes in CTs and inside omg2
wS.omg2(msh.node.CT.onG) = k2*tw(msh.newNodes(msh.node.CT.onG));
%   nodes in CTs and on gam
%   Note that, there may be nodes on Gam but not in NCTs 2. However, we
%   don't care about them because we only consider the NCTs which don't
%   contain them.

%-------------------------------------------------------------------------
wS.ct1 = sparse(msh.nStd,1);
wS.ct1(msh.node.CT.all) = k1*tw(msh.node.CT.all);
%   all nodes in CTs

%-------------------------------------------------------------------------
wS.ct2 = sparse(msh.nStd,1);
wS.ct2(msh.node.CT.all) = k2*tw(msh.newNodes(msh.node.CT.all));
%   all nodes in CTs

end