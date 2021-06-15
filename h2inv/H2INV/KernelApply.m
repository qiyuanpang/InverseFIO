function u = KernelApply(N, f)
if(1)  % This is turned on for the Kernel given implicitly by
       % five_solve.
  global mvcount;
  global Tree;
  global TotalKernelTime;
  mvcount = mvcount + size(f,2);
  Tsta = tic;
  u = five_solve(N, Tree, f);
  Tend = toc(Tsta);
  TotalKernelTime = TotalKernelTime + Tend;
end

if(0) % This is turned on if Kernel is given explicitly by a global
      % matrix Kernel.
  global Kernel;
  u = Kernel * f;
end
