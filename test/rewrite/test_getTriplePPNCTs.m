% %% getGMgPP in main_sys_linda
% % Must run to step 2 because uold=0 at the beginning
% % useNewton = 0;
% 
% K = getKug(tris,CT,uold,msh,pa,cp,typeG);
% [igPP1,jgPP1,vgPP1] = getTriplePPNCTs(NCTs1,msh,pa,K.NC1); % NCTs1
% [igPP2,jgPP2,vgPP2] = getTriplePPNCTs(NCTs2,msh,pa,K.NC2); % NCTs2
% 
% defG = @findG;
% uS1 = uold.omg1;
% uS2 = uold.omg2;
% defH1 = @(x,y,pa) cp.kk1;
% defH2 = @(x,y,pa) cp.kk2;
% [igPP1c,jgPP1c,vgPP1c] = getTriplePPNCTs2(msh,pa,NCTs1,'u',uS1,'gu',defG,'h',defH1);
% [igPP2c,jgPP2c,vgPP2c] = getTriplePPNCTs2(msh,pa,NCTs2,'u',uS2,'gu',defG,'h',defH2);
% 
% % if results are empty double column vectors, then good
% isempty(find(igPP1 ~= igPP1c, 1))
% isempty(find(jgPP1 ~= jgPP1c, 1))
% isempty(find(vgPP1 ~= vgPP1c, 1))
% 
% isempty(find(igPP2 ~= igPP2c, 1))
% isempty(find(jgPP2 ~= jgPP2c, 1))
% isempty(find(vgPP2 ~= vgPP2c, 1))


%% getGMuNewton in main_sys_linda
% Must run to step 2 because uold=0 at the beginning
% useNewton = 1;

K2 = getKwug(tris,CT,uold,wS,msh,pa,typeG);
[iwgPP1,jwgPP1,vwgPP1] = getTriplePPNCTs(NCTs1,msh,pa,K2.NC1); % NCTs1
[iwgPP2,jwgPP2,vwgPP2] = getTriplePPNCTs(NCTs2,msh,pa,K2.NC2); % NCTs2

defG2 = @findG2;
defFw = @findFw;
uS1 = uold.omg1; uS2 = uold.omg2;
wS1 = wS.omg1; wS2 = wS.omg2;

[iwgPP1c,jwgPP1c,vwgPP1c] = getTriplePPNCTs2(msh,pa,NCTs1,...
    'u',uS1,'gu',defG2,'w',wS1,'fw',defFw);
[iwgPP2c,jwgPP2c,vwgPP2c] = getTriplePPNCTs2(msh,pa,NCTs2,...
    'u',uS2,'gu',defG2,'w',wS2,'fw',defFw);

% if results are empty double column vectors, then good
find(iwgPP1 ~= iwgPP1c)
find(jwgPP1 ~= jwgPP1c)
find(vwgPP1 ~= vwgPP1c)

find(iwgPP2 ~= iwgPP2c)
find(jwgPP2 ~= jwgPP2c)
find(vwgPP2 ~= vwgPP2c)


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