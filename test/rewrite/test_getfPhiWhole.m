% %% model: main_sys_linda
% % getPf.m, findDefFtw
% 
% itest = 1;
% triangle = msh.t(1:3,115);
% defF = model.defFtw;
% Ptest = getPf(tris,CT,msh,pa,defF);
% P= Ptest.NC1(3,:);
% Fwold = getLWhole(triangle,itest,msh,pa,P)
% 
% defF = @findDefFtw;
% Fwnew = getfPhiWhole(msh,pa,triangle,itest,'h',defF)
% 
% function valF = findDefFtw(xx,yy,pa)
%     % Define right hand side function f
%     % Input: coordinate of points + indication subdomain
%     % Output: value of phi at points
%     rr2 = (xx-0.5).^2 + (yy-0.5).^2; % r^2
%     valF = -16*rr2; % -16r^2
% end

%% model: main_sys_linda
% test with w,u

typeG = 1; % g(u)=u
itest = 3;
t = 15;
triangle = NCTs1(1:3,t);
Ptest = getPwg(tris,CT,uoldEach,wS,msh,pa,typeG);
P = Ptest.NC1(t,:);
Fwold = getLWhole(triangle,itest,msh,pa,P)

defG = @findG;
defFw = @findFw;
uSn = uoldEach.omg1;
wSn = wS.omg1;
Fwnew = getfPhiWhole(msh,pa,triangle,itest,...
                'u',uSn,'gu',defG,'w',wSn,'fw',defFw)

function val = findG(uu)
    val = uu.^2;
end

function val = findFw(ww)
    val = ww;
end
