% %% model: main_sys_linda
% % getLf of tw
% 
% typeG = 1; % g(u)=u
% NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
% defF = model.defFtw;
% Ptest = getPf(tris,CT,msh,pa,defF);
% % [i1,f1] = getLNCTs(NCTs1,msh,pa,Ptest.NC1); % NCTs1
% [i2,f2] = getLNCTs(NCTs2,msh,pa,Ptest.NC2); % NCTs1
% 
% defH = @findDefH;
% % [in1,fn1] = getfPhiNCTs(msh,pa,NCTs1,'h',defH);
% [in2,fn2] = getfPhiNCTs(msh,pa,NCTs2,'h',defH);
% 
% % if results are empty double column vectors, then good
% % find(i1 ~= in1)
% % find(f1 ~= fn1)
% 
% find(i2 ~= in2)
% find(f2 ~= fn2)
% 
% function valF = findDefH(xx,yy,pa)
%     rr2 = (xx-0.5).^2 + (yy-0.5).^2; % r^2
%     valF = -16*rr2; % -16r^2
% end



%% model: main_sys_linda
% getLfwg with u,w

typeG = 1;
NCTs1=tris.NCTs1; NCTs2=tris.NCTs2;
Ptest = getPwg(tris,CT,uoldEach,wS,msh,pa,typeG);

[i1,f1] = getLNCTs(NCTs1,msh,pa,Ptest.NC1); % NCTs1
[i2,f2] = getLNCTs(NCTs2,msh,pa,Ptest.NC2); % NCTs2

defG = @findG;
defFw = @findFw;

uS1 = uoldEach.omg1;
wS1 = wS.omg1;
[in1,fn1] = getfPhiNCTs(msh,pa,NCTs1,'u',uS1,'gu',defG,'w',wS1,'fw',defFw);

uS2 = uoldEach.omg2;
wS2 = wS.omg2;
[in2,fn2] = getfPhiNCTs(msh,pa,NCTs2,'u',uS2,'gu',defG,'w',wS2,'fw',defFw);


% if results are empty double column vectors, then good
find(i1 ~= in1)
find(f1 ~= fn1)
find(i2 ~= in2)
find(f2 ~= fn2)

function val = findG(uu)
    val = uu.^2;
end

function val = findFw(ww)
    val = ww;
end