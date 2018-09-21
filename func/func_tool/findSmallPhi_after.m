function [trisNew,CTnew] = findSmallPhi_after(msh,pa,phi,tris,CT)
% If there are small cut, change the components to the new ones after
%   removing all small
% Create this function because these of lines of codes are reused frequently
% Status: checked with the old code
% Input: - all triangles: tris
%        - cut triangles and their components: CTs, CT
% Output: new ones


%% Get info
CTs=tris.CTs; NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;


%% initial
trisNew = tris; CTnew = CT;


%% Find small cut
[tri2del,t2Omg1,t2Omg2] = findSmallPhi(msh,CTs,CT,pa,phi);
% If there are small-cut triangles, remove them from CTs!!
if ~isempty(tri2del)
    nCTs = size(CTs,2); % number of OLD cut triangles
    
    % get NEW not-cut triangles
    if ~isempty(t2Omg1)
        NCTs1 = [NCTs1,CTs(:,t2Omg1)]; % add more triangles to NCTs1
        trisNew.NCTs1=NCTs1;
    end
    
    if ~isempty(t2Omg2)
        NCTs2 = [NCTs2,CTs(:,t2Omg2)]; % add more triangles to NCTs2
        trisNew.NCTs2=NCTs2;
    end
    
    % get NEW cut triangles
    CTs = CTs(:,setdiff(1:nCTs,tri2del));
    trisNew.CTs=CTs;
    
    % find again all information
    CTnew = getInfoCTs(CTs,phi,msh,pa);
end

end