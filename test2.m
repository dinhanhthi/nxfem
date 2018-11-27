% find(abs(phi)<pa.tol)
yTop = max(points(2,abs(phi)<pa.tol));
if isempty(yTop)
    tmp(:,:) = iPs(2,:,:); % 2 x nCTs (y coor of all iPs)
    yTop = max(tmp(:));
else
    tmp(:,:) = iPs(2,:,:); % 2 x nCTs (y coor of all iPs)
    yTop = max(yTop,max(tmp(:)));
end

yTop