function [mH11,mH12] = getMatrixH1(omg1NCTs,omg2NCTs,areaNCTs1,...
               areaNCTs2,CTs,iPs,typeCTs,areaChildCTs,areaCTs,idxEachCTs,pa)
% Find the matrix of |e|_H1 in subdomain Omg_i
% Input: not cut triangles, cut triangles
% Output: matrix to be used in E_iA_iE'_i

kk1 = pa.kk1;
kk2 = pa.kk2;

%% ========================================================
% GET TRIPLETS from PHI*PHI terms
% all are column-arrays
% =========================================================
[ip1,jp1,vp1] = getMatrixGijNCTs(omg1NCTs,areaNCTs1,1,pa); % in Omg1
[ip2,jp2,vp2] = getMatrixGijNCTs(omg2NCTs,areaNCTs2,1,pa); % in Omg2
[ipt,jpt,vpt1,vpt2] = getMatrixGijCTs(CTs,iPs,typeCTs,idxEachCTs,...
                                        areaCTs,1,msh,pa); % in CTs



%% ========================================================
% GET TRIPLETS from GRAD*GRAD terms
% all are column-arrays
% =========================================================
[iG1,jG1,vG1] = getMatrixGradGradNCTs(omg1NCTs,kk1,msh); % in Omg1 NCTs
[iG2,jG2,vG2] = getMatrixGradGradNCTs(omg2NCTs,kk2,msh); % in Omg2 NCTs
[iGt,jGt,vGt1,vGt2] = getMatrixGradGradCTs(CTs,areaChildCTs,...
                                    areaCTs,pa,msh); % CTs



%% ========================================================
% BUID MATRIX
% =========================================================
ii1 = [iG1;iGt;ip1;ipt]; 
jj1 = [jG1;jGt;jp1;jpt];
vv1 = [vG1;vGt1;vp1;vpt1];

ii2 = [iG2;iGt;ip2;ipt]; 
jj2 = [jG2;jGt;jp2;jpt];
vv2 = [vG2;vGt2;vp2;vpt2];


mH11 = sparse(ii1,jj1,vv1);
mH12 = sparse(ii2,jj2,vv2);
end