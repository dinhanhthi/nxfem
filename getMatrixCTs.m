function [mCTs1,mCTs2,mCTs3,mCTs4] = getMatrixCTs(CTs,areaChildCTs,...
                            areaCTs,iPs,uNormalsCT,msh,pa,cp)
% Get matrices on cut triangles
% Input: all information about cut triangles
% Output: 4 matrices on cut triangles
%       - A_{ij} for phi_i, phi_j whose support in Omg1 : mCTs1
%       - A_{k(i)k(j)} for phi_i, phi_j whose support in Omg2 : mCTs2
%       - A_{k(i)j} for phi_k(i) in Omg2, phi_j in omg1 : mCTs3
%       - A_{ik(j)} for phi_i in Omg1, phi_k(j) in omg2 : mCTs4

% We consider kk1, kk2 here instead of inside getMatrixGradGradCTs because
%   we are gonna use getMatrixGradGradCTs in another purpose with different
%   kk1 and kk2.

k1=pa.kk1; k2=pa.kk2;

[iGG,jGG,vGG1,vGG2] = getMatrixGradGradCTs(CTs,areaChildCTs,...
                            areaCTs,k1,k2,msh); % CTs (term GradGrad)
[iOG,jOG,vOG1,vOG2,vOG3,vOG4]...
    = getMatrixOnGamCTs(CTs,iPs,uNormalsCT,msh,pa,cp); % CTs (terms on Gam)

% cT1
iCTs1 = [iGG,iOG]; jCTs1 = [jGG,jOG]; vCTs1 = [vGG1,vOG1];
mCTs1 = sparse(iCTs1,jCTs1,vCTs1);
% cT2
iCTs2 = [iGG,iOG]; jCTs2 = [jGG,jOG]; vCTs2 = [vGG2,vOG2];
mCTs2 = sparse(iCTs2,jCTs2,vCTs2);
% cT3
mCTs3 = sparse(iOG,jOG,vOG3);
% cT4
mCTs4 = sparse(iOG,jOG,vOG4);

end