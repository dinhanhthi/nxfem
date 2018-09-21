%% Find the convergence rate for the simple-level-set test cases
% This case: using SUPG and don't apply FMM
% nSeg in [37 57 77 101]
% Related files: main_levelset_simple.m
% Output:
%   - s1: ||phi_h^N-phi^0||_{L2(Omg)}
%   - s2: ||phi^0||_L2(Gam_h^N)
%   - s3: ||phi_h^0-phi^0||_{L2(Omg)}
%   - s4: ||phi^0||_L2(Gam_h^0)


%%
nSeg = [37, 57, 77, 101];

h = [0.066009409376807 0.043592587959105 0.034506407398553 0.025898412629370];

s1 = [0.0057330124 0.0027451846 0.0018080699 0.0011379785];

s2 = [0.0082977 0.0041700 0.0021747 0.0009022];

s3 = [0.0008837343 0.0002934948 0.0001892726 0.0001223398];
s4 = [0.0010321 0.0003997 0.0002475 0.0001297];

% for i=2:4
%     order1 = log(s1(i-1)/s1(i))/log(h(i-1)/h(i));
%     display(order1);
% end
% 
% for i=2:4
%     order2 = log(s2(i-1)/s2(i))/log(h(i-1)/h(i));
%     display(order2);
% end

for i=2:4
    order3 = log(s3(i-1)/s3(i))/log(h(i-1)/h(i));
    display(order3);
end

% for i=2:4
%     order4 = log(s4(i-1)/s4(i))/log(h(i-1)/h(i));
%     display(order4);
% end

plot(log(h),log(s1),'b',log(h),log(s2),'r',log(h),log(s3),'g',log(h),log(s4),'k');

legend('s1','s2','s3','s4')

% n1 = polyfit(log(h),log(s1),1)
% n2 = polyfit(log(h),log(s2),1)
% n3 = polyfit(log(h),log(s3),1)
% n4 = polyfit(log(h),log(s4),1)