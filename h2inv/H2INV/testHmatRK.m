% TESTHMATRK is a test subroutine.  It generates a random low-rank
% Kernel.  This kernel is then compressed using H^2 compression.  The
% two kernels are compared using a random test vector.
%
% NOTE: To use this subroutine, KernelApply.m should be changed
% accordingly.  See KernelApply.m for more details.

if(1) 
  global Kernel;
  N = 32;
  ntotal = N*N;
  GRD = reshape(1:ntotal,N,N);

  Kernel = zeros(ntotal, ntotal);

  nLevel = 3;
  display('Constructing Sample List first');
  tic

    SampleList = ConstructSampleList( N, nLevel );

  toc
  display(['Constructing Random Kernel upto Level ' int2str(nLevel)]);
  tic

    nRealRank = 3;
    nWidth = N / 4;
    for iLevel = 1 : nLevel
      if( iLevel > 1 )
	nWidth = nWidth / 2;
      end
      for iSource = 1 : SampleList{iLevel}.nSource
	SourceBlock = SampleList{iLevel}.Source(:, iSource);
	Source = [(SourceBlock(1)-1)*nWidth+1, ...
	  (SourceBlock(1))*nWidth; ...
	  (SourceBlock(2)-1)*nWidth+1, ...
	  (SourceBlock(2))*nWidth];
	for iTarget = 1 : SampleList{iLevel}.nTarget
	  TargetBlock = SampleList{iLevel}.SourceTarget(:, iTarget, ...
	    iSource );
	  Target = [(TargetBlock(1)-1)*nWidth+1, ...
	    (TargetBlock(1))*nWidth; ...
	    (TargetBlock(2)-1)*nWidth+1, ...
	    (TargetBlock(2))*nWidth];
	  SourceInd = reshape( GRD(Source(1,1):Source(1,2), ...
	    Source(2,1):Source(2,2)), nWidth * nWidth, 1 );
	  TargetInd = reshape( GRD(Target(1,1):Target(1,2), ...
	    Target(2,1):Target(2,2)), nWidth * nWidth, 1 );
	  Kernel(SourceInd, TargetInd) = ...
	    randn(nWidth*nWidth,nRealRank) * ...
	    randn(nRealRank,nWidth*nWidth) + ...
	    1d-10 * randn(nWidth*nWidth, nWidth * nWidth );
	end
      end
    end

    Kernel = (Kernel + transpose(Kernel))/2;
  toc
end

if(1)
  display(['Constructing H matrix upto Level ' int2str(nLevel)]);
  tic

    epsilon = 1e-6;
    svdfact = 1;
    [DiagVec, DiagBlock, nDiagBlock, SampleList, RKTree] = ...
      ConstructHMatrix( N, nLevel, epsilon, svdfact );
  toc

end


if(1)
  display(['Verifying the H-matrix tree' int2str(nLevel)]);
  tic

  u = rand(N, N);
  vec = reshape( Kernel * reshape(u, ntotal, 1), N, N );
  vecRK = RKTreeApplyH2( SampleList, RKTree, u, N, nLevel+1, 1 );


  toc
  norm( vec - vecRK, 'fro' ) ./ norm( vec, 'fro' )
end
