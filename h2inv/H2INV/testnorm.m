function [TotalTime, TotalKernelTime, mvcount, TotalBytes, L2abserrpower, L2relerrpower] = ...
  testnorm(N, nLevel)
% TESTNORM is a test subroutine for the peeling algorithm for H^2
% matrix.  The Kernel is an opearator of the form (-Delta + V)^{-1}. 
%
% INPUT parameter:
%   N            :   grid size per dimension.
%   nLevel       :   number of low rank level.
%
% OUTPUT parameter:
%   TotalTime    :   Total computational time.
%   mvcount      :   Total number of matrix-vector multiplication.
%   TotalBytes   :   Total memory cost.
%   L2abserrpower:   Absolute L^2 error of the Green's function
%                    estimated by power method.
%   L2relerrpower:   Relative L^2 error of the Green's function
%                    estimated by power method.
%
% AVAILABLE FUNCTIONS:
%  (*)  Test error for diagonal blocks
%  (*)  Estimate 2-norm of the Green's function using random vectors.
%  (*)  Estimate 2-norm of the Green's function using power method.
%  (*)  Estimate the memory cost.
if(1)

  h = (1/N);
  V = 1 + rand(N,N);
  % V = 1e-3*ones(N,N);

  global Tree;
  disp('Factorization');
  Tree = five_setup(V,N,h);
  disp('Extraction');
  DiagExact = five_extract(N,Tree);
end


if(1)
  svdfact = 1;
  epsilon = 1e-6;

  global mvcount TotalTime TotalKernelTime;
  mvcount = 0;
  TotalKernelTime = 0;
  global DiagVec  DiagBlock  nDiagBlock  SampleList  RKTree;
  [DiagVec, DiagBlock, nDiagBlock, SampleList, RKTree] = ...
    ConstructHMatrix( N, nLevel, epsilon, svdfact );
  disp(['Number of matrix-vector applications is ' num2str(mvcount)])
  fprintf('\n')

  disp(['relative L^1 error of the diagonal elements']);
  diagrelerrL1 = sum(sum(abs(DiagVec - DiagExact))) ./ ...
    sum(sum(abs(DiagExact)))

  disp(['absolute Max error of the diagonal elements']);
  diagabserrMax = max(max(abs(DiagVec - DiagExact)))

  disp(['relative Max error of the diagonal elements']);
  diagrelerrMax = max(max(abs(DiagVec - DiagExact))) ./ ...
    max(max(abs(DiagExact)))

end

% Estimate 2-norm using random method
if(0)
  nSample = 100;
  u = rand(N,N,nSample);
  v1 = WholeTreeApply( ...
    SampleList, RKTree, DiagBlock, u, N, nSample, nLevel+1, nDiagBlock);
  v2 = KernelApply( N, reshape(u, N*N, nSample) );
  disp(['Estimated L^2 error using ' num2str(nSample) ' tests']);
  L2err = max( sqrt(sum((v1-v2).^2,1)) ./ sqrt(sum(v1.^2,1)) )
end



% Estimate 2-norm using power method. The criterion is to measure the
% eigenvalue difference.
if(1)
  nSample = 1;
  nPowerMax = 500;
  u = rand(N,N,nSample);
  v = reshape(u, N*N, nSample);
  Eold = 1;
  Enew = 0;
  for iter = 1 : nPowerMax
    if( mod(iter, 20) == 0 )
      fprintf('%5i th power iteration for (G-G_approx)\n\n', iter);
    end
    w = reshape(v, N, N, nSample);
    v1 = WholeTreeApply( ...
      SampleList, RKTree, DiagBlock, w, N, nSample, nLevel+1, nDiagBlock);
    v2 = KernelApply( N, reshape(w, N*N, nSample) );
    vdiff = reshape(v1-v2, N*N, nSample);
    Enew = abs(sum(v.* vdiff, 1));
    if( abs(Enew - Eold)./ Enew < 1e-3 )
      fprintf('%5i iterations to reach abs(Enew - Eold)./ Enew < 1e-3 \n\n', iter);
      break;
    end
    v = vdiff;
    Eold = Enew;
    for j = 1 : nSample
      v(:,j) = v(:,j) / norm(v(:,j),2);
    end 
  end
  L2abserrpower = Enew
  disp(['Estimated absolute L^2 error using power method']);
  L2abserrpower
end

% Estimate the L^2 norm of G matrix to calculate relative L^2 error
if(1)
  nSample = 1;
  nPowerMax = 500;
  u = rand(N,N,nSample);
  v = reshape(u, N*N, nSample);
  Eold = 1;
  Enew = 0;
  for iter = 1 : nPowerMax
    if( mod(iter, 20) == 0 )
      fprintf('%5i th power iteration for (G-G_approx)\n\n', iter);
    end
    w = v;
    v1 = KernelApply( N, reshape(w, N*N, nSample) );
    Enew = abs(sum(v1.*v,1));
    if( abs(Enew - Eold)./ Enew < 1e-3 )
      fprintf('%5i iterations to reach abs(Enew - Eold)./ Enew < 1e-3 \n\n', iter);
      break;
    end
    Eold = Enew;
    v = v1;
    for j = 1 : nSample
      v(:,j) = v(:,j) / norm(v(:,j),2);
    end 
  end
  L2G = Enew;
  disp(['Estimated relative L^2 error using power method']);
  L2relerrpower = L2abserrpower / L2G
end

%%
% Estimate memory cost for both the low rank tree and the diagonal and
% near-off-diagonal blocks for H^2 matrix. We only count the real
% numbers. The indices are not counted.  
if(1)
  TotalBytes = 0;
  S = RKTree{nLevel}.RKBlockH2U;
  infoS = whos('S');
  TotalBytes = TotalBytes + infoS.bytes;
  S = RKTree{nLevel}.RKBlockH2M;
  infoS = whos('S');
  TotalBytes = TotalBytes + infoS.bytes;
  for iLevel = 1 : nLevel-1
    S = RKTree{iLevel}.RKBlockH2T;
    infoS = whos('S');
    TotalBytes = TotalBytes + infoS.bytes;
    S = RKTree{iLevel}.RKBlockH2M;
    infoS = whos('S');
    TotalBytes = TotalBytes + infoS.bytes;
  end
  for iDiagBlock = 1 : nDiagBlock
    S = DiagBlock{iDiagBlock}.Mat;
    infoS = whos('S');
    TotalBytes = TotalBytes + infoS.bytes;
  end
  fprintf('H2: Total Memory Cost for RKTree and DiagBlock (real numbers only) is %6.3f MB\n', ...
    TotalBytes/1024^2);
end


%%
% Estimate memory cost for both the low rank tree and the diagonal and
% near-off-diagonal blocks for Uniform H^1 matrix.
if(1)
  TotalBytesH1 = 0;
  for iLevel = 1 : nLevel
    S = RKTree{iLevel}.RKBlockUniU;
    infoS = whos('S');
    TotalBytesH1 = TotalBytesH1 + infoS.bytes;
    S = RKTree{iLevel}.RKBlockUniM;
    infoS = whos('S');
    TotalBytesH1 = TotalBytesH1 + infoS.bytes;
  end
  for iDiagBlock = 1 : nDiagBlock
    S = DiagBlock{iDiagBlock}.Mat;
    infoS = whos('S');
    TotalBytesH1 = TotalBytesH1 + infoS.bytes;
  end
  fprintf('Uniform H^1: Total Memory Cost for RKTree and DiagBlock (real numbers only) is %6.3f MB\n', ...
    TotalBytesH1/1024^2);
end
