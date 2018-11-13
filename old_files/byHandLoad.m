ct=2;
d1=3;
d2=1;
d3=2;
kk = 1;
xuoi=1; % cái này ip2 tiep ngay sau dinh tam giác chu ko phai là xuôi cua iP1-iP2
tugiac=1;
if xuoi
    ip1 = iPs(:,1,ct);
    ip2 = iPs(:,2,ct);
else
    ip2 = iPs(:,1,ct);
    ip1 = iPs(:,2,ct);
end
if tugiac
    tgb = rhsTGct(ct,CTs);
    tgn = rhsPart(d1,d2,d3,ip1,ip2,ct,CTs,kk);
    kq = tgb-tgn;
else
   kq = rhsPart(d1,d2,d3,ip1,ip2,ct,CTs,kk);
end

t=4;
val=rhsTGt(t)

kq

% kq1 = kq
% kq2 = kq
% kq3 = kq

function val = rhsPart(d1,d2,d3,ip1,ip2,ct,CTs,kk)
global points
Stg = getAreaTri(points(:,CTs(d1,ct)),...
                points(:,CTs(d2,ct)),points(:,CTs(d3,ct)));

[~,y1] = getCoorRef(ip1,points(:,CTs(d1,ct)),...
                points(:,CTs(d2,ct)),points(:,CTs(d3,ct)));
[x2,~] = getCoorRef(ip2,points(:,CTs(d1,ct)),...
                points(:,CTs(d2,ct)),points(:,CTs(d3,ct)));
switch kk
    case 1
       val = -(2/3)*Stg*(3*y1*x2-y1*x2^2-y1^2*x2); 
    case 2
        val = -(2/3)*Stg*x2^2*y1;
    case 3
        val = -(2/3)*Stg*y1^2*x2;
end
end

function val=rhsTGt(t)
global triangles points
Stg = getAreaTri(points(:,triangles(1,t)),...
                points(:,triangles(2,t)),points(:,triangles(3,t)));
val = -2/3*Stg;
end

function val=rhsTGct(ct,CTs)
global points
Stg = getAreaTri(points(:,CTs(1,ct)),...
                points(:,CTs(2,ct)),points(:,CTs(3,ct)));
val = -2/3*Stg;
end