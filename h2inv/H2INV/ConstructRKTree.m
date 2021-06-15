function RKTree = ConstructRKTree( N, nLevel, SampleList )
% ********************************************************
% INPUT PARAMETERS
% ********************************************************
% N         :   Domain size is N*N
% nLevel    :   The number of low rank levels
% SampleList:   Cell structure, containing the information of the
%               decomposd domain.
% ********************************************************
% OUTPUT PARAMETERS
% ********************************************************
% RKTree    :   Cell structure, containing the low rank matrices and
%               related information.  For more details:
%               nWidth    : Size of each block at iLevel.
%               nBlock    : Number of blocks at iLevel.
%               nRKBlock  : Number of low rank matrices at iLevel. 
%                           nRKBlock = nTarget*nSource/2
%               RKSource  : 2*2*nRKBlock array. The ABSOLUTE INDEX of
%                           the Source in 2D domain. For instance, if
%                           N=16, RKTree{1}.RKSource(:,:,1) =
%                           [1 8; 1 8]; 
%                           The first line is the starting and ending
%                           index in row (target) direction, and the
%                           second line is the starting and ending index
%                           in column (source) direction.
%               RKTarget  : 2*2*nRKBlock array. Similar to RKSource, the
%                           ABSOLUTE INDEX of the Target in 2D domain.
%               RKSourceInd: nRKBlock*1 array. The index for the source
%                            of the low rank matrix in SampleList.Source.
%               RKTargetInd: nRKBlock*1 array. The index for the target
%                            of the low rank matrix in SampleList.Source.
%               RKBlockUniU: nSource*1 cell. U matrices after
%                            uniformization.
%               RKBlockUniM: nRKBlock*1 cell. M matrices (a.k.a Sigma
%                            matrices) that is compatible with U after
%                            uniformization.
%               RKH2S     : The diagonal of S matrices in the
%                            uniform H^1 representation.
%               RKBlockH2U:  nSource*1 cell. U matrices after H^2
%                            compression.
%               RKBlockH2M:  nRKBlock*1 cell. M matrices that is
%                            compatible with U after H^2 compression.
%               RKBlockH2T:  nChild * nSource cell.  Transfer matrix T
%                            that connects the H2U at current level (or
%                            more exactly, H2UChild) and H2U of each of
%                            its child.
%               RKBlockH2UChild:  nChild * nSource cell.  H2U matrix at
%                                 current level restricted to each of
%                                 its child.  This is used in H^2
%                                 compression.
%               RKBlockU   : nRKBlock*1 cell. U matrices before
%                            uniformization. Will be freed after the
%                            uniformization step.
%               RKBlockV   : nRKBlock*1 cell. V matrices before
%                            uniformization. Will be freed after the
%                            uniformization step.
%               RKBlockM   : nRKBlock*1 cell. M matrices before
%                            uniformization. Will be freed after the
%                            uniformization step.
%               RKBlockSDIRECT   : nRKBlock*1 cell. Singular values in
%                                  the SVD step in peeling-off process.
%                                  This block is combined with RKBlockU
%                                  during uniformization. Will be freed
%                                  after uniformization.
%               RKBlockSINDIRECT : nRKBlock*1 cell. Singular values in
%                                  the SVD step in peeling-off process.
%                                  This block is combined with RKBlockV
%                                  during uniformization. Will be freed
%                                  after uniformization.
%               RKRankU          : The rank of U matrices for each
%                                  low rank matrix. For debug use.
%               RKRankV          : The rank of V matrices for each
%                                  low rank matrix. For debug use.
%               RKRankUniU       : The rank of uniformized U matrices.
%                                  For debug use.
%               RKRankH2U        : The rank of H^2 U matrices. For debug use.

DIRECT = 1;
INDIRECT = 0;

RKTree = cell(nLevel,1);
nWidthTop = N / 4;

for iLevel = 1 : nLevel
  if( iLevel == 1 )
    nWidth = nWidthTop;
    nBlockCur = 4;
  else
    nWidth = nWidth / 2;
    nBlockCur = nBlockCur * 2;
  end
  nSourceCur = SampleList{iLevel}.nSource;
  nTargetCur = SampleList{iLevel}.nTarget;
  nRKBlockCur = nSourceCur * nTargetCur / 2;
  
  nChild = SampleList{iLevel}.nChild;

  % Each element of RKSource is a 2 by 2 matrix recording the position
  % (in 2-d domain) of the source; Similarly for RKTarget
  
  RKTree{iLevel} = struct('nWidth', nWidth, ...
    'nBlock', nBlockCur, 'nRKBlock', nRKBlockCur, ...
    'RKSource', [], 'RKTarget', [], ...
    'RKSourceInd', [], 'RKTargetInd', [], ...
    'RKBlockU', [], 'RKBlockV', [], ...
    'RKBlockUniU', [], ...
    'RKH2S', [], ...
    'RKBlockUniM', [], ...
    'RKBlockH2U', [], ...
    'RKBlockH2T', [], ...
    'RKBlockH2M', [], ...
    'RKBlockH2UChild', [], ...
    'RKBlockM', [], ...
    'RKBlockSDIRECT', [], ...
    'RKBlockSINDIRECT', [], ...
    'RKRankU', [], ...
    'RKRankV', [], ...
    'RKRankUniU', [], ...
    'RKRankH2U', [] );

  RKTree{iLevel}.RKSource = zeros( 2, 2, RKTree{iLevel}.nRKBlock );
  RKTree{iLevel}.RKTarget = zeros( 2, 2, RKTree{iLevel}.nRKBlock );
  RKTree{iLevel}.RKSourceInd = zeros( RKTree{iLevel}.nRKBlock, 1 );
  RKTree{iLevel}.RKTargetInd = zeros( RKTree{iLevel}.nRKBlock, 1 );
  RKTree{iLevel}.RKBlockU = cell(nRKBlockCur, 1);
  RKTree{iLevel}.RKBlockV = cell(nRKBlockCur, 1);
  RKTree{iLevel}.RKBlockUniU = cell(nSourceCur, 1);
  RKTree{iLevel}.RKH2S = cell(nSourceCur, 1);
  RKTree{iLevel}.RKBlockM = cell(nRKBlockCur, 1);
  RKTree{iLevel}.RKBlockUniM = cell(nRKBlockCur, 1);

  RKTree{iLevel}.RKBlockH2U = cell(nSourceCur, 1);
  RKTree{iLevel}.RKBlockH2T = cell(nChild, nSourceCur);
  RKTree{iLevel}.RKBlockH2M = cell(nRKBlockCur, 1);
  RKTree{iLevel}.RKBlockH2UChild = cell(nChild, nSourceCur);

  RKTree{iLevel}.RKBlockSDIRECT = cell(nRKBlockCur, 1);
  RKTree{iLevel}.RKBlockSINDIRECT = cell(nRKBlockCur, 1);
 
  RKTree{iLevel}.RKRankU = zeros(nRKBlockCur,1);
  RKTree{iLevel}.RKRankV = zeros(nRKBlockCur,1);
  RKTree{iLevel}.RKRankUniU = zeros(nSourceCur,1);
  RKTree{iLevel}.RKRankH2U = zeros(nSourceCur,1);



  for iSource = 1 : SampleList{iLevel}.nSource
    for iTarget = 1 : SampleList{iLevel}.nTarget
      if( SampleList{iLevel}.SymFlag( iTarget, iSource ) == ...
	  DIRECT )
	Source = SampleList{iLevel}.Source(:, iSource);
	Target = SampleList{iLevel}.SourceTarget(:, iTarget, iSource );
	RKTree{iLevel}.RKSourceInd(SampleList{iLevel}.RKIndex(...
	  iTarget, iSource )) = iSource;
	RKTree{iLevel}.RKTargetInd(SampleList{iLevel}.RKIndex(...
	  iTarget, iSource )) = ...
	  SampleList{iLevel}.SourceTargetInd(iTarget, iSource);
	% The following 2*2 matrix is organized as
	% [ y_0 y_1 ]
	% [ x_0 x_1 ]
	RKTree{iLevel}.RKSource( :, :, SampleList{iLevel}.RKIndex(...
	  iTarget, iSource ) ) = ...
	  [ 1+(Source(1)-1)*RKTree{iLevel}.nWidth, ...
	    Source(1)*RKTree{iLevel}.nWidth; 
	    1+(Source(2)-1)*RKTree{iLevel}.nWidth, ...
	    Source(2)*RKTree{iLevel}.nWidth ];
	RKTree{iLevel}.RKTarget( :, :, SampleList{iLevel}.RKIndex(...
	  iTarget, iSource ) ) = ...
	  [ 1+(Target(1)-1)*RKTree{iLevel}.nWidth, ...
	    Target(1)*RKTree{iLevel}.nWidth; 
	    1+(Target(2)-1)*RKTree{iLevel}.nWidth, ...
	    Target(2)*RKTree{iLevel}.nWidth;
	  ];
      end
    end
  end
end
