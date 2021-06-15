function SampleList = ConstructSampleList(N, nLevel)
% ********************************************************
% INPUT PARAMETERS
% ********************************************************
% N         :   Domain size is N*N
% nLevel    :   The number of low rank levels
% ********************************************************
% OUTPUT PARAMETERS
% ********************************************************
% SampleList:   Cell structure of length nLevel, containing the
%               information of the decomposd domain. For more details:
%               nSource:    Number of sources at iLevel.
%               nTarget:    Number of targets associated to EACH source
%                           iLevel. For the current case,
%                           nTarget = 7     if iLevel == 1,
%                           nTarget = 27    if iLevel >  1.
%               nParallel:  Number of sources that can be combined
%                           together at iLevel.
%               Source:     2*nSource array. The BLOCK POSITION of each
%                           source at iLevel. For example, at first
%                           level the domain is partitioned into 4*4
%                           blocks, then the block position ranges from
%                           [1;1] to [4;4].  NOTE: Source is NOT ordered
%                           in the natural order starting from Level 3.
%               Target:     2*nTarget array. The RELATIVE BLOCK POSITION
%                           of each target. The real target block
%                           position is the source position plus the
%                           relative target position.
%               RKIndex:    nTarget*nSource array. All the low rank
%                           matrices are saved in a one dimensional
%                           array in RKTree. RKIndex(iTarget, iSource)
%                           returns this index corresponds to the
%                           (iTarget, iSource) pair. Note that G(Target,
%                           Source) and G(Source, Target) corresponds to
%                           the same block.
%               SourceTarget: 2*nTarget*nSource array. The BLOCK
%                           POSITION of the a Target corresponding to a Source.
%               SourceTargetInd:  nTarget*nSource array. The index of a
%                                 Target corresponding to a Source in
%                                 the array SampleList{iLevel}.Source. 
%               Symflag:    nTarget*nSource array. Whether the low rank
%                           matrix (Target, Source) is a DIRECT or an
%                           INDIRECT matrix.
%               nChild:     Scallar.  Number of childs for each level.
%                           This is currently set to 4.
%               Child:      nChild * nSource array.  The indices of the
%                           childs for a source at current level.  NOTE:
%                           This index corresponds to the index in
%                           SampleList{nextlevel}.Source instead of
%                           natural order.


SampleList = cell(nLevel, 1);
DIRECT = 1;
INDIRECT = 0;

for iLevel = 1:nLevel
  if( iLevel == 1 )
    nBlockCur = 4; 
    nSourceCur = 16;
    nTargetCur = 7;
    nParallel = 1; 
    % Number of sources that we can work together; the sources are
    % stored in the way that 1 .. nParallel, nParallel+1
    % .. 2*nParallel, ... can be combined together.

    nChild = 4;

    SampleList{iLevel} = struct('nSource', nSourceCur, 'nTarget', ...
                                nTargetCur, 'nParallel', nParallel, ...
                                'Source', [], 'Target', [], ...
				'RKIndex', [], 'SourceTarget', [], ...
                                'SourceTargetInd', [], ...
				'SymFlag', [], 'nChild', nChild, ...
                                'Child',[] );

    SampleList{iLevel}.Child = zeros( nChild, SampleList{iLevel}.nSource );

    SampleList{iLevel}.Source = zeros( 2, SampleList{iLevel}.nSource);
    SampleList{iLevel}.Target = zeros( 2, SampleList{iLevel}.nTarget);


    for iSource = 1:nSourceCur
      SampleList{iLevel}.Source(:, iSource) = [mod(iSource-1, nBlockCur)+1, ...
                    ceil(iSource/nBlockCur)];
    end % iSource 
    

  

    %FORTRAN convention is used for labelling.
                       
    for iTarget = 1:nBlockCur
      SampleList{iLevel}.Target(:, iTarget) = [iTarget, 2];
    end
    for iTarget = nBlockCur+1:nTargetCur
      SampleList{iLevel}.Target(:, iTarget) = [2, iTarget - nBlockCur - 2];
    end % iTarget

  else
    nBlockCur = nBlockCur * 2;
    nSourceCur = nSourceCur * 4; 
    nTargetCur = 27; 
    if (iLevel > 2) 
      nParallel = nParallel * 4;
      nParaBlock = nParaBlock * 2;
    else
      nBlockPrev = nBlockCur;
      nParaBlock = 1;
    end

    nChild = 4;

    SampleList{iLevel} = struct('nSource', nSourceCur, 'nTarget', ...
                                nTargetCur, 'nParallel', nParallel, ...
                                'Source', [], 'Target', [], ...
				'RKIndex', [], 'SourceTarget', [], ...
                                'SourceTargetInd', [], ...
				'SymFlag', [], 'nChild', nChild, ...
                                'Child',[] );

    SampleList{iLevel}.Child = zeros( nChild, SampleList{iLevel}.nSource );
    
    SampleList{iLevel}.Source = zeros( 2, SampleList{iLevel}.nSource);
    SampleList{iLevel}.Target = zeros( 2, SampleList{iLevel}.nTarget, ...
                                       4);

    for iSource = 1 : round(nSourceCur/nParallel)
      for iParallel = 1 : nParallel
        SampleList{iLevel}.Source(:, (iSource-1) * nParallel + ...
            iParallel) = [mod(iSource-1, nBlockPrev)+1 + ...
                    (mod(iParallel-1, nParaBlock)) * nBlockPrev, ...
                    ceil(iSource/nBlockPrev) + ...
                    (ceil(iParallel/nParaBlock) - 1) * nBlockPrev];
      end
    end % iSource
   

    % Target is the same for level 2 and above. Since we are saving the
    % (Target, Source) pair, there is actually no need to store
    % individually the target arrays.
    iTarget = 0;
    for jShift = [2 3 -2]
      for iShift = [-2:3]
        iTarget = iTarget + 1;
        SampleList{iLevel}.Target(:, iTarget, 1) = [iShift, jShift];
      end
    end
    for jShift = [-1:1]
      for iShift = [2 3 -2]
        iTarget = iTarget + 1;
        SampleList{iLevel}.Target(:, iTarget, 1) = [iShift, jShift];
      end
    end
    
    iTarget = 0;
    for jShift = [2 3 -2]
      for iShift = [-3:2]
        iTarget = iTarget + 1;
        SampleList{iLevel}.Target(:, iTarget, 2) = [iShift, jShift];
      end
    end
    for jShift = [-1:1]
      for iShift = [2 -3 -2]
        iTarget = iTarget + 1;
        SampleList{iLevel}.Target(:, iTarget, 2) = [iShift, jShift];
      end
    end

    iTarget = 0;
    for jShift = [2 -3 -2]
      for iShift = [-2:3]
        iTarget = iTarget + 1;
        SampleList{iLevel}.Target(:, iTarget, 3) = [iShift, jShift];
      end
    end
    for jShift = [-1:1]
      for iShift = [2 3 -2]
        iTarget = iTarget + 1;
        SampleList{iLevel}.Target(:, iTarget, 3) = [iShift, jShift];
      end
    end

    iTarget = 0;
    for jShift = [2 -3 -2]
      for iShift = [-3:2]
        iTarget = iTarget + 1;
        SampleList{iLevel}.Target(:, iTarget, 4) = [iShift, jShift];
      end
    end
    for jShift = [-1:1]
      for iShift = [2 -3 -2]
        iTarget = iTarget + 1;
        SampleList{iLevel}.Target(:, iTarget, 4) = [iShift, jShift];
      end
    end
    
  end % if( iLevel == 1 )
 
  % Construct the mapping between (Target, Source)  
  cIndex = 0;
  SampleList{iLevel}.SourceTarget = ...
    zeros( 2, nTargetCur, nSourceCur );
  SampleList{iLevel}.SourceTargetInd = ...
    zeros( nTargetCur, nSourceCur );
  if( iLevel == 1 )
    for iSource = 1 : nSourceCur
      Source = SampleList{iLevel}.Source(:, iSource);
      for iTarget = 1 : nTargetCur
	Target = Source + SampleList{iLevel}.Target(:, iTarget);
	Target = mod( Target-1, nBlockCur ) + 1;
	SampleList{iLevel}.SourceTarget( :, iTarget, iSource ) = Target;
	[tempL, SampleList{iLevel}.SourceTargetInd(iTarget, iSource)] = ...
          ismember(Target', SampleList{iLevel}.Source', 'rows');
      end
    end
  else
    for iSource = 1 : nSourceCur
      Source = SampleList{iLevel}.Source(:, iSource);
      RowKind = mod(Source(1)-1, 2)+1;
      ColKind = mod(Source(2)-1, 2)+1;
      TargetKind = RowKind + (ColKind-1)*2;
      for iTarget = 1 : nTargetCur
	Target = Source + SampleList{iLevel}.Target(:, iTarget, ...
	  TargetKind);
	Target = mod( Target-1, nBlockCur ) + 1;
	SampleList{iLevel}.SourceTarget( :, iTarget, iSource ) = Target;
	[tempL, SampleList{iLevel}.SourceTargetInd(iTarget, iSource)] = ...
          ismember(Target', SampleList{iLevel}.Source', 'rows');
      end
    end
  end

  % Construct the RK indicies
 
  SampleList{iLevel}.RKIndex = zeros( nTargetCur, nSourceCur );
  % SymFlag == DIRECT cooresponds to lower-diagonal blocks
  % SymFlag == INDIRECT cooresponds to upper-diagonal blocks, because it
  % is the tranposes of these blocks that we want.
   
  SampleList{iLevel}.SymFlag = zeros( nTargetCur, nSourceCur );
  
  for iSource = 1 : nSourceCur
    Source = SampleList{iLevel}.Source(:, iSource);
    indSource = Source(1) + (Source(2)-1) * nBlockCur;
    for iTarget = 1 : nTargetCur
      Target = SampleList{iLevel}.SourceTarget(:, iTarget, ...
	iSource);
      indTarget = Target(1) + (Target(2)-1) * nBlockCur;
      if( indSource < indTarget )
	cIndex = cIndex + 1;
	SampleList{iLevel}.RKIndex( iTarget, iSource ) = cIndex;	
	SampleList{iLevel}.SymFlag( iTarget, iSource ) = DIRECT;
      end
    end
  end
  if( cIndex ~= nSourceCur * nTargetCur / 2 )
    error('Something wrong with the number of boxes');
  end

   
  for iSource = 1 : nSourceCur
    Source = SampleList{iLevel}.Source(:, iSource);
    indSource = Source(1) + (Source(2)-1) * nBlockCur;
    for iTarget = 1 : nTargetCur
      Target = SampleList{iLevel}.SourceTarget(:, iTarget, ...
	iSource );
      indTarget = Target(1) + (Target(2)-1) * nBlockCur;
      if( indSource > indTarget )
	indtmp1 = find( SampleList{iLevel}.Source(1, :)  ...
	  == Target(1) );
	indtmp2 = find( SampleList{iLevel}.Source(2, indtmp1) ...
	  == Target(2) );
	iSourceReverse = indtmp1(indtmp2);
	indtmp1 = find( SampleList{iLevel}.SourceTarget(1, :, ...
	  iSourceReverse ) == Source(1) );
	indtmp2 = find( SampleList{iLevel}.SourceTarget(2, indtmp1, ...
	  iSourceReverse ) == Source(2) );
	iTargetReverse = indtmp1(indtmp2);
	SampleList{iLevel}.RKIndex( iTarget, iSource ) = ...
	  SampleList{iLevel}.RKIndex( iTargetReverse, iSourceReverse);
	SampleList{iLevel}.SymFlag( iTarget, iSource ) = INDIRECT;
      end
    end
  end

end % iLevel

% Assign Child arrays
% NOTE: Starting from Level 3, the domain is no longer indexed using the
% natural order.  Therefore when assigning the child indices for each
% domain at current level, a search step is mandatory to obtain a
% consistent index number!
for iLevel = 1 : nLevel-1
  for iSource = 1:SampleList{iLevel}.nSource
    Row = SampleList{iLevel}.Source(1,iSource);
    Col = SampleList{iLevel}.Source(2,iSource);
    RowChild = 2*Row-1;  ColChild = 2*Col-1;
    [Tmp1, Tmp2] = ismember([RowChild ColChild], SampleList{iLevel+1}.Source', 'rows');
    SampleList{iLevel}.Child(1, iSource) = Tmp2;
    RowChild = 2*Row;  ColChild = 2*Col-1;
    [Tmp1, Tmp2] = ismember([RowChild ColChild], SampleList{iLevel+1}.Source', 'rows');
    SampleList{iLevel}.Child(2, iSource) = Tmp2;
    RowChild = 2*Row-1;  ColChild = 2*Col;
    [Tmp1, Tmp2] = ismember([RowChild ColChild], SampleList{iLevel+1}.Source', 'rows');
    SampleList{iLevel}.Child(3, iSource) = Tmp2;
    RowChild = 2*Row;  ColChild = 2*Col;
    [Tmp1, Tmp2] = ismember([RowChild ColChild], SampleList{iLevel+1}.Source', 'rows');
    SampleList{iLevel}.Child(4, iSource) = Tmp2;
  end % iSource 
end


