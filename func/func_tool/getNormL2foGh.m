function val = getNormL2foGh(msh,pa,CTs,iPs,defF)
% Find ||f(x,y)||_L^2(Gam_h) where f is a function handle
% Note: first used in finding find norm of exact phi on Gam_h, like in (7.48) p.224
% Arnold Book
% Input: - defF(x,y,pa) : function handle
%        - CTs, iPs: cut triangles and inter section points
% Output: value of norm

nCTs = size(CTs,2); % number of cut triangles
points = msh.p;

dim=1; deg=pa.degP1D; % quadrature's info
[wt,pt] = getGaussQuad(dim,deg); 
nwt = size(wt,2); % number of Gaussian points

val = 0;
for t = 1:nCTs
    pointA = iPs(:,1,t); % 1st intersection point
    pointB = iPs(:,2,t); % 2nd intersection point

    v1 = points(:,CTs(1,t)); % vertx 1
    v2 = points(:,CTs(2,t)); % vertx 2
    v3 = points(:,CTs(3,t)); % vertx 3

    lenAB = sqrt((pointB(1)-pointA(1))^2+(pointB(2)-pointA(2))^2);

    [xHa,yHa] = getCoorRef(pointA,v1,v2,v3);
    [xHb,yHb] = getCoorRef(pointB,v1,v2,v3);

    tmp = 0;
    for k = 1:nwt
        xHaT = xHa + (xHb-xHa)*((1+pt(k))/2); % xHa(t), t=(1+xi)/2
        yHaT = yHa + (yHb-yHa)*((1+pt(k))/2); % yHa(t), t=(1+xi)/2
        [xk,yk] = getCoorSTD([xHaT,yHaT],v1,v2,v3);
        tmp = tmp + wt(k)*defF(xk,yk,pa)*defF(xk,yk,pa);
    end
    
    val = val + 0.5*lenAB*tmp;
end

val = sqrt(val);

end