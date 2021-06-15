% MAINTEST is a test subroutine that performs multiple test comparison
% and exports results in LATEX format.  This subroutines depends on
% testnorm.m.
 
nTest = 1;

N = zeros(nTest,1);
nLevel = zeros(nTest,1);
TotalTime = zeros(nTest,1);
TotalKernelTime = zeros(nTest,1);
mvcount = zeros(nTest,1);
TotalBytes = zeros(nTest,1);
L2abserrpower = zeros(nTest, 1);
L2relerrpower = zeros(nTest, 1);

% N = [32 64 64];
% nLevel = [3 3 4];

N = [32];
nLevel = [3];

for i = 1 : nTest
  [TotalTime(i), TotalKernelTime(i), mvcount(i), TotalBytes(i), L2abserrpower(i), ...
    L2relerrpower(i)] = testnorm(N(i), nLevel(i)) 
end

table1 = [];
table1 = [table1 sprintf(' \\sqrt{N} & L & Total Time (sec) & Construction Time (sec) & matvec number & Memory (MB)\\\\ \n')];
table1 = [table1 sprintf(' \\hline\n')];
for i = 1 : nTest
  table1 = [table1 sprintf('%4d & %2d & %6.2f & %6.2f %6d & %6.2f\\\\ \n', ...
    N(i), nLevel(i), TotalTime(i), TotalTime(i) - TotalKernelTime(i), mvcount(i), TotalBytes(i)/1024^2)];
end
table1 = [table1 sprintf(' \\hline\n')];

table2 = [];
table2 = [table2 sprintf(' \\sqrt{N} & L & Absolute L^2 error & Relative L^2 error\\\\ \n')];
table2 = [table2 sprintf(' \\hline\n')];
for i = 1 : nTest
  table2 = [table2 sprintf('%4d & %2d & %10.2e & %10.2e \\\\ \n', ...
    N(i), nLevel(i), L2abserrpower(i), L2relerrpower(i))];
end
table2 = [table2 sprintf(' \\hline\n')];


table1
table2
