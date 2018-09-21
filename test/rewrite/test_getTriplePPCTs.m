%% getGMgPP in main_sys_linda
K = getKug(tris,CT,uold,msh,pa,cp,typeG);
[igPPc,jgPPc,vgPPc1,vgPPc2] = getTriplePPCTs(CTs,CT,msh,pa,K); % CTs

defG = @findG;
uS1 = uold.ct1;
uS2 = uold.ct2;
defH1 = @(x,y,pa) cp.kk1;
defH2 = @(x,y,pa) cp.kk2;
[igPPc1,jgPPc1,vgPPc1t] = getTriplePPCTs2(msh,pa,CTs,CT,1,'u',uS1,'gu',defG,'h',defH1);
[igPPc2,jgPPc2,vgPPc2t] = getTriplePPCTs2(msh,pa,CTs,CT,2,'u',uS2,'gu',defG,'h',defH2);


% if results are empty double column vectors, then good
find(igPPc ~= igPPc1)
find(jgPPc ~= jgPPc1)
find(vgPPc1 ~= vgPPc1t)

find(igPPc ~= igPPc2)
find(jgPPc ~= jgPPc2)
find(vgPPc2 ~= vgPPc2t)



%% getGMuNewton in main_sys_linda
K2 = getKwug(tris,CT,uold,wS,msh,pa,typeG);
[iwgPPc,jwgPPc,vwgPPc1,vwgPPc2] = getTriplePPCTs(CTs,CT,msh,pa,K2); % CTs

defG2 = @findG2;
defFw = @findFw;
uS1 = uold.ct1; uS2 = uold.ct2;
wS1 = wS.ct1; wS2 = wS.ct2;

[iwgPPc1,jwgPPc1,vwgPPc1t] = getTriplePPCTs2(msh,pa,CTs,CT,1,...
    'u',uS1,'gu',defG2,'w',wS1,'fw',defFw);
[iwgPPc2,jwgPPc2,vwgPPc2t] = getTriplePPCTs2(msh,pa,CTs,CT,2,...
    'u',uS2,'gu',defG2,'w',wS2,'fw',defFw);

% if results are empty double column vectors, then good
find(iwgPPc ~= iwgPPc1)
find(jwgPPc ~= jwgPPc1)
find(vwgPPc1 ~= vwgPPc1t)

find(iwgPPc ~= iwgPPc2)
find(jwgPPc ~= jwgPPc2)
find(vwgPPc2 ~= vwgPPc2t)


%% function handles
function val = findG(uu)
    val = uu.^2;
end
function val = findG2(uu)
    val = 2.*uu;
end
function val = findFw(ww)
    val = ww;
end