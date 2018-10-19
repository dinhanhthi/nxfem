%% Find the convergence rate for the simple-level-set test cases
% This case: using both SUPG and FMM
% nSeg in [37 57 77 101]
% Related files: main_levelset_simple.m
% Output:
%   - s1: ||phi_h^N-phi^0||_{L2(Omg)}
%   - s2: ||phi^0||_L2(Gam_h^N)
%   - s3: ||phi_h^0-phi^0||_{L2(Omg)}
%   - s4: ||phi^0||_L2(Gam_h^0)

%
nSeg = [37, 57, 77, 101];
h = [0.0660094 0.04359258 0.03450640739 0.0258984126];
s1 = [0.0633557017 0.0495139390 0.0502913544 0.0511053598];
s2 = [0.0100336 0.0055870 0.0034501 0.0019325];

% s3 = [0.0008837343 0.0002934948 0.0001892726 0.0001223398];

for i=2:4
    order1 = log(s1(i-1)/s1(i))/log(h(i-1)/h(i));
    display(order1);
end

for i=2:4
    order2 = log(s2(i-1)/s2(i))/log(h(i-1)/h(i));
    display(order2);
end

plot(log(h),log(s1),'b',log(h),log(s2),'r');

legend('s1','s2')

% n1 = polyfit(log(h),log(s1),1)
% n2 = polyfit(log(h),log(s2),1)
% n3 = polyfit(log(h),log(s3),1)