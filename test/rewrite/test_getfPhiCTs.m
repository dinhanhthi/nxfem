% %% model: main_sys_linda
% % test with h(x,y)
% 
% defF = model.defFtw;
% Ptest = getPf(tris,CT,msh,pa,defF);
% [ic,fc1,fc2] = getLCTs(CTs,CT,msh,pa,Ptest); % CTs
% 
% defH1 = @(x,y,pa) defF(x,y,pa,1);
% [in1,fn1] = getfPhiCTs(msh,pa,CTs,CT,1,'h',defH1);
% defH2 = @(x,y,pa) defF(x,y,pa,2);
% [in2,fn2] = getfPhiCTs(msh,pa,CTs,CT,2,'h',defH2);
% 
% % if results are empty double column vectors, then good
% find(ic ~= in1)
% find(fc1 ~= fn1)
% find(ic ~= in2)
% find(fc2 ~= fn2)

%% model: main_sys_linda
% test with u,w

typeG = 1; % cf. defG.m
Ptest = getPwg(tris,CT,uoldEach,wS,msh,pa,typeG);
[ic,fc1,fc2] = getLCTs(CTs,CT,msh,pa,Ptest); % CTs

defGu = @findGu;
defFw = @findFw;
uS1 = uoldEach.ct1; wS1 = wS.ct1;
uS2 = uoldEach.ct2; wS2 = wS.ct2;
[in1,fn1] = getfPhiCTs(msh,pa,CTs,CT,1,'u',uS1,'gu',defGu,'w',wS1,'fw',defFw);
[in2,fn2] = getfPhiCTs(msh,pa,CTs,CT,2,'u',uS2,'gu',defGu,'w',wS2,'fw',defFw);

% if results are empty double column vectors, then good
find(ic ~= in1)
find(fc1 ~= fn1)
find(ic ~= in2)
find(fc2 ~= fn2)

function val = findGu(uu)
    val = uu.^2;
end

function val = findFw(ww)
    val = ww;
end
