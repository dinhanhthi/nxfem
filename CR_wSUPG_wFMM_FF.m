%% Find the convergence rate for the simple-level-set test cases
% This case: - using both SUPG and FMM
%            - LIMIT number of using FMM
%            - mesh taken from freefem++ (for comparison with sol obtained
%            from it), h follows exp 2
% Related files: main_levelset_simple.m
% Output:
%   - s1: ||phi_h^N-phi^0||_{L2(Omg)}
%   - s2: ||phi^0||_L2(Gam_h^N)
%   - s3: ||phi_h^0-phi^0||_{L2(Omg)}
%   - s4: ||phi^0||_L2(Gam_h^0)

h = [0.0977536973 0.0481075631 0.0279129707];
s1 = [0.0761112661 0.0652214143 0.0667357419];
s2 = [0.0271405356 0.0093166425 0.0019122083];


for i=2:3
    order1 = log(s1(i-1)/s1(i))/log(h(i-1)/h(i));
    display(order1);
end

for i=2:3
    order2 = log(s2(i-1)/s2(i))/log(h(i-1)/h(i));
    display(order2);
end

plot(log(h),log(s1),'b',log(h),log(s2),'r');
legend('s1','s2')