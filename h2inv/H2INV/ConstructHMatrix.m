function [DiagVec, DiagBlock, nDiagBlock, SampleList, RKTree] = ...
  ConstructHMatrix(N, nLevel, epsilon, svdfact)
% ********************************************************
% INPUT PARAMETERS
% ********************************************************
% N         :   Domain size is N*N
% nLevel    :   The number of low rank levels
% epsilon   :   the SVD thresholding criterion at the bottom level
% svdfact   :   the growth factor of SVD thresholding criterion between
%               each levels:
%               SVDCut(iLevel) = SVDCut(iLevel-1) * svdfact;
%
% ********************************************************
% OUTPUT PARAMETERS
% ********************************************************
% DiagVec   :   N*N array. The diagonal elements of Green's function.
% DiagBlock :   Cell structure, containing the part of Green's function
%               that is not represented by low rank matrics. This part
%               includes the diagonal blocks as well as off-diagonal
%               blocks.
% SampleList:   Cell structure, containing the information of the
%               decomposd domain.
% RKTree    :   Cell structure, containing the low rank matrices and
%               related information. 
% nDiagBlock:   The length of DiagBlock.
%
% ********************************************************
% OTHER TERMINOLOGIES
% ********************************************************
% G         :   short for Green's function.
% Source    :   The column direction of Green's function. This comes
%               from the matrix vector multiplication, where the vector
%               acts as sources and the column direction multiplies with
%               the sources.
% Target    :   The row direction of Green's function. This comes
%               from the matrix vector multiplication, after which the
%               resulting vector is the contribution of the sources to
%               the targets.
% DIRECT    :   If G(Target, Source) is in the lower diagonal matrix,
%               G(Target, Source) is called DIRECT. During the
%               peeling-off process, we put test sources in Source
%               region, and extract information from Target region.
%               Since the matrix G is symmetric, the lower diagonal
%               matrices provide the U part in the SVD decomposition of
%               G(Target,Source), and the upper diagonal matrices
%               provide the V part in the SVD decomposition.
% INDIRECT  :   If G(Target, Source) is in the upper diagonal matrix,
%               G(Target, Source) is called INDIRECT.               

global TotalTime;
TotalTime = 0;

DIRECT = 1;
INDIRECT = 0;

nWidthTop = N / 4;
% nWidthBottom = 32;
% nLevel = round( log2( nWidthTop / nWidthBottom ) + 1);

maxRank = zeros(nLevel, 1);
maxRank(1) = 20;
maxRank(2:end) = 10;

nSample = maxRank + 10;
SVDCut = zeros(nLevel, 1);
SVDCut(1) = epsilon;
for iLevel = 2 : nLevel
  SVDCut(iLevel) = SVDCut(iLevel-1) * svdfact;
end
SVDCut

SampleList = ConstructSampleList( N, nLevel );
RKTree = ConstructRKTree(N, nLevel, SampleList);

display('Constructing RKTree...');
if(1)
  % Peeling-off for level 1 to nLevel.
  for iLevel = 1 : nLevel
    display( iLevel )
    Tstart = tic;

    nWidth = RKTree{iLevel}.nWidth;
    nGroup = SampleList{iLevel}.nSource / ...
      SampleList{iLevel}.nParallel;
    TempR = randn( N, N, nSample(iLevel) );
    TempQ = cell(RKTree{iLevel}.nRKBlock, 1);
    % TempR: random sources.
    % TempQ: saving target information after matrix vector
    % multiplication for each low-rank matrices.

    for i = 1 : RKTree{iLevel}.nRKBlock
      TempQ{i} = zeros( nWidth, nWidth, nSample(iLevel) );
    end

    % Applying the kernel to the random matrices
    for iGroup = 1 : nGroup

      TempMatApply = zeros( N, N, nSample(iLevel) );
      for iParallel = 1 : SampleList{iLevel}.nParallel
	iSource = (iGroup-1) * SampleList{iLevel}.nParallel + ...
	  iParallel;
	Source = SampleList{iLevel}.Source( :, iSource );
	TempMatApply( 1+(Source(1)-1)*nWidth : ...
	  Source(1)*nWidth, ...
	  1+(Source(2)-1)*nWidth : ...
	  Source(2)*nWidth, : ) = ...
	  TempR( 1+(Source(1)-1)*nWidth : ...
	  Source(1)*nWidth, ...
	  1+(Source(2)-1)*nWidth : ...
	  Source(2)*nWidth, : );
      end
      TempMatRes = KernelApply( N, reshape(...
	TempMatApply, N * N, nSample(iLevel) ) );
      TempMatRes = reshape( TempMatRes, N, N, nSample(iLevel) ); 
      TempMatRKRes = RKTreeApplyH2( SampleList, RKTree, TempMatApply, N, iLevel, nSample(iLevel) );
      TempMatRes = TempMatRes - TempMatRKRes;

      % Extracting the blocks for each target and perform SVD
      for iParallel = 1 : SampleList{iLevel}.nParallel
	iSource = (iGroup-1) * SampleList{iLevel}.nParallel + ...
	  iParallel;
	for iTarget = 1 : SampleList{iLevel}.nTarget
	  iRKBlock = SampleList{iLevel}.RKIndex(iTarget, iSource );
	  Target = SampleList{iLevel}.SourceTarget( :, ...
	    iTarget, iSource );
	  TempT = TempMatRes(...
	    1+(Target(1)-1)*nWidth : ...
	    Target(1)*nWidth,  ...
	    1+(Target(2)-1)*nWidth : ...
	    Target(2)*nWidth, : );
	  [TempU, TempS, TempW ] = svd( ...
	    reshape( TempT, nWidth * nWidth, nSample(iLevel) ), 'econ' );
	  nRankU = min( max(find( abs(diag(TempS)) ./ max(abs(diag(TempS))) > SVDCut(iLevel) )), ...
	    maxRank(iLevel) );
	  if( isempty(nRankU) )
	    nRankU = 0;
	  end
	  if( nRankU == maxRank(iLevel)  )
	    disp('max rank reached')
	  end
	  % DIRECT matrices contribute to U part and Sigma part, INDIRECT
	  % matrices contribute to V part.
	  if( SampleList{iLevel}.SymFlag(iTarget, iSource ) == DIRECT )
	    RKTree{iLevel}.RKBlockU{ iRKBlock } = TempU(:, 1:nRankU );
	    RKTree{iLevel}.RKBlockSDIRECT{ iRKBlock } = ...
	      (TempS(1:nRankU,1:nRankU));
	    TempQ{iRKBlock} = TempT;
	    RKTree{iLevel}.RKRankU(iRKBlock) = nRankU;
	  else
	    RKTree{iLevel}.RKBlockV{ iRKBlock } = TempU(:, 1:nRankU );
	    RKTree{iLevel}.RKBlockSINDIRECT{ iRKBlock } = ...
	      (TempS(1:nRankU,1:nRankU));
	    RKTree{iLevel}.RKRankV(iRKBlock) = nRankU;
	  end
	  clear TempT TempU TempS TempW;
	end % iTarget
      end % iParallel	

      clear TempMatApply TempMatRes TempMatRKRes;

    end % iGroup

    % Calculating Sigma and form low rank approximation
    for iRKBlock = 1 : RKTree{iLevel}.nRKBlock
      Source = RKTree{iLevel}.RKSource(:, :, iRKBlock);
      Target = RKTree{iLevel}.RKTarget(:, :, iRKBlock);
      TempR1 = reshape( ...
	TempR( Source(1,1):Source(1,2), Source(2,1):Source(2,2), : ), ...
	nWidth * nWidth, nSample(iLevel) );
      TempR2 = reshape( ...
	TempR( Target(1,1):Target(1,2), Target(2,1):Target(2,2), : ), ...
	nWidth * nWidth, nSample(iLevel) );

      RKTree{iLevel}.RKBlockM{iRKBlock} = pinv( transpose(TempR2) * ...
	RKTree{iLevel}.RKBlockU{iRKBlock} ) * ...
	transpose(TempR2) * reshape( ...
	TempQ{iRKBlock}, nWidth * nWidth, nSample(iLevel) ) * ...
	pinv( transpose(RKTree{iLevel}.RKBlockV{iRKBlock}) * ...
	TempR1 );

      clear TempR1 TempR2;
    end % iRKBlock

    clear TempR TempQ;

    % Uniformization
    for iSource = 1 : SampleList{iLevel}.nSource  
      UAll = [];
      for iTarget = 1 : SampleList{iLevel}.nTarget
	iRKBlock = SampleList{iLevel}.RKIndex(iTarget, iSource );
	% Note that the following line may look very strange at a first
	% glance. Remember Source is the column direction and Target is
	% the row direction for the matrix. SymFlag is DIRECT if this is
	% a lower diagonal block matrix. However, G_{ij} need to be
	% represented as
	%    G_{i,j} = U_{i} M_{ij} U_{j}^T 
	%            = u_{ij} m_{ij} v_{ij}^T.
	% Therefore we need to compress V part when matrix is DIRECT,
	% and compress U part when matrix is INDIRECT.
	%
	% Weighted uniformization is performed here for performance
	% reason.
	if( SampleList{iLevel}.SymFlag(iTarget, iSource) == ...
	    DIRECT )
	  UAll = [UAll RKTree{iLevel}.RKBlockV{iRKBlock} * ...
	    RKTree{iLevel}.RKBlockSINDIRECT{iRKBlock}];
	else
	  UAll = [UAll RKTree{iLevel}.RKBlockU{iRKBlock} * ...
	    RKTree{iLevel}.RKBlockSDIRECT{iRKBlock}];
	end
      end %iTarget

      % Perform SVD to get U after uniformization. In principle this
      % should be done using random SVD. Since this is not very
      % computational demanding, we use regular SVD decomposition here
      % first.
      [TempU, TempS, TempV] = svd(UAll, 'econ');

      nRankU = ...
	min( max(find( abs(diag(TempS)) ./ max(abs(diag(TempS))) > SVDCut(iLevel) )), ...
	4*maxRank(iLevel) );
      % The number 4 is just a empirical upper bound after
      % uniformization, usually the factor is much less than 4.
      if( isempty(nRankU) )
	nRankU = 0;
      end
      if( nRankU == maxRank(iLevel) *4  )
	disp('max rank reached by uniformization')
      end
      RKTree{iLevel}.RKBlockUniU{iSource} = TempU(:, 1:nRankU);
      RKTree{iLevel}.RKRankUniU(iSource) = nRankU; 
      % Save S matrix for H^2 compression use.
      RKTree{iLevel}.RKH2S{iSource} = diag(TempS(1:nRankU, 1:nRankU));
    end %iSource

    nGroupWidth = round(sqrt(SampleList{iLevel}.nSource / ...
      SampleList{iLevel}.nParallel));
    nParallelWidth = round(sqrt(SampleList{iLevel}.nParallel));
    nParallel = SampleList{iLevel}.nParallel;


    % Recalculate Sigma matrices (M matrices here) after uniformization.
    for iSource = 1 : SampleList{iLevel}.nSource
      for iTarget = 1 : SampleList{iLevel}.nTarget
	iRKBlock = SampleList{iLevel}.RKIndex(iTarget, iSource );
	if( SampleList{iLevel}.SymFlag(iTarget, iSource) == ...
	    DIRECT )

	  TargetInd = RKTree{iLevel}.RKTargetInd(iRKBlock);

	  % For DIRECT matrices, 
	  % G(Target,Source) = U(Target) * M(Target,Source) * U(Source)'
	  RKTree{iLevel}.RKBlockUniM{iRKBlock} = ...
	    transpose(RKTree{iLevel}.RKBlockUniU{TargetInd}) * ...
	    RKTree{iLevel}.RKBlockU{iRKBlock} * ...
	    RKTree{iLevel}.RKBlockM{iRKBlock} * ...
	    transpose(RKTree{iLevel}.RKBlockV{iRKBlock}) * ...;
	  RKTree{iLevel}.RKBlockUniU{iSource};
	end
      end  % iTarget
    end % iSource

    % Construct H2 matrix

    if( iLevel > 1 )
      for iSourcePre = 1 : SampleList{iLevel-1}.nSource  
	for iChild = 1 : SampleList{iLevel-1}.nChild
	  iSourceCur = SampleList{iLevel-1}.Child(iChild,iSourcePre);
	  if(1)
	     % UAll = [RKTree{iLevel}.RKBlockUniU{iSourceCur}  ...
	     % RKTree{iLevel-1}.RKBlockH2UChild{iChild,iSourcePre}];

	     % size(RKTree{iLevel}.RKBlockUniU{iSourceCur})
	     % size(diag(RKTree{iLevel}.RKH2S{iSourceCur}))
	     % size(RKTree{iLevel-1}.RKBlockH2UChild{iChild,iSourcePre})
	     % size(diag(RKTree{iLevel-1}.RKH2S{iSourcePre}))
	    
	     UAll = [RKTree{iLevel}.RKBlockUniU{iSourceCur} * ...
	     diag(RKTree{iLevel}.RKH2S{iSourceCur}) ...
	     RKTree{iLevel-1}.RKBlockH2UChild{iChild,iSourcePre} * ...
	     diag(RKTree{iLevel-1}.RKH2S{iSourcePre})];

	    % Perform SVD to construct H^2 representation from uniform H^1
	    % representation. In principle this should be done using
	    % random SVD. Since this is not very computational demanding,
	    % we use regular SVD decomposition here first.
	    [TempU, TempS, TempV] = svd(UAll, 'econ');

	    nRankU = ...
	      min( max(find( abs(diag(TempS)) ./ max(abs(diag(TempS))) > SVDCut(iLevel) )), ...
	      4*maxRank(iLevel) );

	    RKTree{iLevel}.RKH2S{iSourceCur} = diag(TempS(1:nRankU, 1:nRankU));
	    RKTree{iLevel}.RKRankH2U(iSourceCur) = nRankU; 

	    % size(RKTree{iLevel}.RKBlockUniU{iSourceCur})
	    % size(UAll)
	    % nRankU
	    % pause

	    % The number 4 is just a empirical upper bound after
	    % uniformization, usually the factor is much less than 4.
	    if( nRankU == maxRank(iLevel) *4  )
	      disp('max rank reached by H^2 construction')
	    end

	    TempU = TempU(:, 1:nRankU);
	  end
	  
	  RKTree{iLevel}.RKBlockH2U{iSourceCur} = TempU;
	  RKTree{iLevel-1}.RKBlockH2T{iChild,iSourcePre} = ...
	    transpose(TempU)* RKTree{iLevel-1}.RKBlockH2UChild{iChild,iSourcePre};
	  % TempU*transpose(TempU)*RKTree{iLevel}.RKBlockUniU{iSourceCur} - RKTree{iLevel}.RKBlockUniU{iSourceCur}
	end %iChild
      end %iSourcePre


      % Recalculate Sigma matrices (M matrices here) after H^2
      % construction.
      for iRKBlock = 1 : RKTree{iLevel}.nRKBlock
	SourceInd = RKTree{iLevel}.RKSourceInd(iRKBlock);
	TargetInd = RKTree{iLevel}.RKTargetInd(iRKBlock);
	RKTree{iLevel}.RKBlockH2M{iRKBlock} = ...
	  transpose(RKTree{iLevel}.RKBlockH2U{TargetInd}) * ...
	  RKTree{iLevel}.RKBlockUniU{TargetInd} * ...
	  RKTree{iLevel}.RKBlockUniM{iRKBlock} * ...
	  transpose(RKTree{iLevel}.RKBlockUniU{SourceInd}) * ...
	  RKTree{iLevel}.RKBlockH2U{SourceInd};
      end % iRKBlock
    end %if (iLevel > 1)

    if (iLevel == 1)
      for iSource = 1 : SampleList{iLevel}.nSource
	RKTree{iLevel}.RKBlockH2U{iSource} = ...
	  RKTree{iLevel}.RKBlockUniU{iSource};
      end
      for iRKBlock = 1 : RKTree{iLevel}.nRKBlock
	RKTree{iLevel}.RKBlockH2M{iRKBlock} = ...
	  RKTree{iLevel}.RKBlockUniM{iRKBlock};
      end 
    end % if (iLevel==1) 

    if( iLevel < nLevel )
      nWidthCur = RKTree{iLevel}.nWidth;
      nWidthChild = RKTree{iLevel+1}.nWidth;
      tmpG = reshape((1:nWidthCur*nWidthCur), nWidthCur, nWidthCur);
      for iSource = 1 : SampleList{iLevel}.nSource
	IND = reshape(tmpG(1:nWidthChild, 1:nWidthChild),[],1);
	RKTree{iLevel}.RKBlockH2UChild{1, iSource} = ...
	  RKTree{iLevel}.RKBlockH2U{iSource}(IND,:);
	IND = reshape(tmpG(nWidthChild+1:2*nWidthChild, 1:nWidthChild),[],1);
	RKTree{iLevel}.RKBlockH2UChild{2, iSource} = ...
	  RKTree{iLevel}.RKBlockH2U{iSource}(IND,:);
	IND = reshape(tmpG(1:nWidthChild, nWidthChild+1:2*nWidthChild),[],1);
	RKTree{iLevel}.RKBlockH2UChild{3, iSource} = ...
	  RKTree{iLevel}.RKBlockH2U{iSource}(IND,:);
	IND = reshape(tmpG(nWidthChild+1:2*nWidthChild, nWidthChild+1:2*nWidthChild),[],1);
	RKTree{iLevel}.RKBlockH2UChild{4, iSource} = ...
	  RKTree{iLevel}.RKBlockH2U{iSource}(IND,:);
      end
    end % iLevel < nLevel


    % Clear memory
    RKTree{iLevel}.RKBlockU = [];
    RKTree{iLevel}.RKBlockV = [];
    RKTree{iLevel}.RKBlockM = [];
    RKTree{iLevel}.RKBlockSDIRECT = [];
    RKTree{iLevel}.RKBlockSINDIRECT = [];

    % Clear memory
    % RKTree{iLevel}.RKBlockUniU = [];
    % RKTree{iLevel}.RKBlockUniM = [];


    Tend = toc(Tstart);
    disp(['Time elasped is ' num2str(Tend) ' s'])
    TotalTime = TotalTime + Tend;
  end % iLevel
end

% Extract diagonal and nearest off diagonal blocks 

fprintf('\n');
disp('Extracting diagonal and nearest off diagonal blocks...');
Tstart = tic;

% If only diagonal blocks are to be extracted, nWidth = nWidth * 2.
% However, we are also extracting the nearest off-diagonal blocks, so
% the increasing factor has to be larger than 2. Optimally this number
% would be 3, but in practice 4 is easier to implement. For simplity
% reason, we quadruple nWidth to sample blocks in parallel. 

nWidthPre = nWidth;
% nWidthPre becomes the width of nLevel-th block size.
nWidth = nWidth * 4;
nBlockCur = round( N / nWidth );

if( 1 ) 
  % The row dimension variables : Target, i, k, or Y
  % The column dimension variables: Source, j, l, or X

  nDiagBlock = 4*4*nBlockCur*nBlockCur;
  DiagBlock = cell(nDiagBlock, 1);
  % Here we use a funny way to save diagonal and nearest off-diagonal
  % blocks as a series of dense blocks, instead of using a general
  % sparse matrix. This is the searching process is too expensive
  % to construct such a sparse matrix using the MATLAB functions.
  for i = 1 : nDiagBlock
    DiagBlock{i} = struct('indSource', [], ...
      'indTarget', [], ...
      'Mat', [] );
  end
  % Diagonal elements of Green's function. 
  DiagVec = zeros(N,N);

  cDiagBlock = 0;
  for j = 1 : 4
    for i = 1 : 4
      for k = 1 : nBlockCur
	for l = 1 : nBlockCur
	  iSource = (k-1)*nWidth + (i-1)*nWidthPre + 1;
	  jSource = (l-1)*nWidth + (j-1)*nWidthPre + 1;
	  iSourceArr = (iSource:iSource+nWidthPre-1)';
	  jSourceArr = (jSource:jSource+nWidthPre-1);
	  indSource = repmat( (jSourceArr-1)*N, length(iSourceArr), 1 ) + ...
	    repmat( iSourceArr, 1, length(jSourceArr) );
	  indSource = reshape(indSource, [], 1);
	  iSourceBlk = (k-1)*4 + i;
	  jSourceBlk = (l-1)*4 + j;
	  iTarget1 = (iSourceBlk-2) * nWidthPre + 1;
	  iTarget2 = (iSourceBlk+1) * nWidthPre;
	  jTarget1 = (jSourceBlk-2) * nWidthPre + 1;
	  jTarget2 = (jSourceBlk+1) * nWidthPre;
	  iTargetArr = (iTarget1:iTarget2)';
	  jTargetArr = jTarget1:jTarget2;
	  iTargetArr = mod( iTargetArr - 1, N ) + 1;
	  jTargetArr = mod( jTargetArr - 1, N ) + 1;
	  indTarget = repmat( (jTargetArr-1)*N, length(iTargetArr), 1 ) + ...
	    repmat( iTargetArr, 1, length(jTargetArr) );
	  indTarget = reshape(indTarget, [], 1);
	  cDiagBlock = cDiagBlock + 1;
	  DiagBlock{cDiagBlock}.indSource = indSource;
	  DiagBlock{cDiagBlock}.indTarget = indTarget;
	end
      end
    end
  end

  % nVecMax is the number of right hand sides that are bundled to solve
  % at one time. This number is determined by the memory size.
  nVecMax = 50;

  TempMatApply = zeros( N, N, nVecMax );
  TempMatRes = zeros(N, N, nVecMax);
  TempMatRKRes = zeros(N, N, nVecMax);
  cVec = 0;
  arrSourceY = zeros(nVecMax);
  arrSourceX = zeros(nVecMax);
  for iSourceX = 1 : nWidth
    for iSourceY = 1 : nWidth

      cVec = cVec + 1;
      arrSourceY(cVec) = iSourceY;
      arrSourceX(cVec) = iSourceX;
      % Parallel treatment
      SourceYInd = iSourceY : nWidth : N;
      SourceXInd = iSourceX : nWidth : N;
      TempMatApply(SourceYInd, SourceXInd, cVec) = 1;
      if( (cVec < nVecMax) && ((iSourceX<nWidth) || (iSourceY<nWidth)) )
	% Wait and solve nVecMax vectors together.
	continue;
      end

      TempMatRes = KernelApply( N, reshape(...
	TempMatApply, N * N, nVecMax) );
      TempMatRes = reshape( TempMatRes, N, N, nVecMax ); 
      TempMatRKRes = RKTreeApplyH2( SampleList, RKTree, TempMatApply, N, nLevel+1, nVecMax);
      TempMatRes = TempMatRes - TempMatRKRes;
      TempMatRes = reshape( TempMatRes, N*N, nVecMax);

      % cVec might not be the same as nVecMax for the last few vectors.
      for iVec = 1 : cVec
	YSourceBlk = floor((arrSourceY(iVec)-1)/nWidthPre) + 1;
	XSourceBlk = floor((arrSourceX(iVec)-1)/nWidthPre) + 1;
	indSourceBlk = (XSourceBlk-1)*4 + YSourceBlk;
	YBlock = mod( arrSourceY(iVec)-1, nWidthPre ) + 1;
	XBlock = mod( arrSourceX(iVec)-1, nWidthPre ) + 1;
	% indBlock is the index inside the block.
	indBlock = (XBlock-1)*nWidthPre + YBlock;
	% Record the diagonal and nearest off-diagonal blocks for all
	% the parallel treated blocks.
	for l = 1 : nBlockCur
	  for k = 1 : nBlockCur
	    cDiagBlock = nBlockCur*nBlockCur*(indSourceBlk-1) + ...
	      (k-1)*nBlockCur + l;
	    DiagBlock{cDiagBlock}.Mat(:, indBlock) = ...
	      TempMatRes(DiagBlock{cDiagBlock}.indTarget, iVec);
	    YVec = arrSourceY(iVec) + (k-1) * nWidth;
	    XVec = arrSourceX(iVec) + (l-1) * nWidth;
	    DiagVec(YVec,XVec) = ...
	      TempMatRes(YVec + (XVec-1)*N, iVec);
	  end
	end
      end

      cVec = 0;
      TempMatApply = zeros( N, N, nVecMax );
      TempMatRes = zeros(N, N, nVecMax);
      TempMatRKRes = zeros(N, N, nVecMax);

    end
  end

  clear TempMatApply TempMatRes TempMatRKRes

end

Tend = toc(Tstart);
disp(['Time elasped is ' num2str(Tend) ' s'])
TotalTime = TotalTime + Tend;

fprintf('\n');
disp(['Total time cost for H-matrix construction is ' ...
  num2str(TotalTime) ' s']);
fprintf('\n');
% whos
