%% Find the convergence rate for the simple-level-set test cases
% This case: - don't use either SUPG nor FMM
%            - mesh taken from freefem++ (for comparison with sol obtained
%            from it), h follows exp 2
% Related files: main_levelset_simple.m
% Output:
%   - s1: ||phi_h^N-phi^0||_{L2(Omg)}
%   - s2: ||phi^0||_L2(Gam_h^N)


% h = [0.1964185496 0.0977536973 0.0481075631 0.0279129707];

h = [0.1382446030	0.0692133200	0.0406824766	0.0210005271];

s1 = [0.0034055608 0.0009189743 0.0003133805 0.0000655538];
s2 = [0.0049372115 0.0009990582 0.0003254385 0.0000751663];


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