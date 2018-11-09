tic; a=0;
A = rand(12000, 4400);
B = rand(12000, 4400);
a=toc-a;
disp(a)


tic; a=0;
C = A'.*B';
% a=toc-a;
disp(toc-a)