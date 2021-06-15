function v = RKTreeApplyH2( SampleList, RKTree, u, N, nLevelCur, nSample )
% RKTreeApplyH2 calculates the low rank matrix-vector multiplication
% where the matrix is H^2.
% ********************************************************
% INPUT PARAMETERS
% ********************************************************
% SampleList:   Cell structure, containing the information of the
%               decomposd domain.
% RKTree    :   Cell structure, containing the low rank matrices and
%               related information. 
% u         :   N * N * nSample array (3D). Vectors to be multiplied.
%               Interpreted as sources.
% N         :   Domain size is N*N
% nLevelCur :   The current peeling-off level. This number equals to (the
%               length of RKTree) + 1.
% nSample   :   Number of right hand sides.
% ********************************************************
% OUTPUT PARAMETERS
% ********************************************************
% v         :   N * N * nSample array (3D). Vectors after matrix-vector
%               multiplcation.  Interpreted as targets.
% 
DIRECT = 1;
INDIRECT = 0;

v = zeros( N, N, nSample );
if( nLevelCur == 1 )
  return;
end


TempT = cell(SampleList{nLevelCur-1}.nSource, nLevelCur-1);
TempS = cell(SampleList{nLevelCur-1}.nSource, nLevelCur-1);

iLevel = nLevelCur - 1;
nWidth = RKTree{iLevel}.nWidth;

% iRow is used because here it does not represent the column (source) direction
% of the matrix, 
for iRow = 1:SampleList{iLevel}.nSource
  Row = SampleList{iLevel}.Source( :, iRow);
  TempT{iRow, iLevel} = transpose(RKTree{iLevel}.RKBlockH2U{iRow}) * ...
    reshape( u( 1+(Row(1)-1)*nWidth : Row(1)*nWidth, ...
    1+(Row(2)-1)*nWidth : Row(2)*nWidth, : ), ...
    nWidth*nWidth, nSample );
end

for iLevel = nLevelCur -2 : -1 : 1
  for iSource = 1:SampleList{iLevel}.nSource
    TempT{iSource, iLevel} = zeros( ...
        size(RKTree{iLevel}.RKBlockH2T{1,iSource},2), nSample);
    for iChild = 1:SampleList{iLevel}.nChild
      
      % % Identify the dimension mismatch
      % if ((nLevelCur == 4) )
	% iLevel
	% iSource
	% iChild
	% size(TempT{iSource, iLevel})
        % size(RKTree{iLevel}.RKBlockH2T{1, iSource})
        % size(TempT{SampleList{iLevel}.Child(iChild, iSource), iLevel+1})
      % end
      
      TempT{iSource, iLevel} = TempT{iSource, iLevel} + ...
          transpose(RKTree{iLevel}.RKBlockH2T{iChild, iSource}) * ...
          TempT{SampleList{iLevel}.Child(iChild, iSource), iLevel+1};
    end
  end
end


for iLevel = 1 : nLevelCur - 1

  for iSource = 1 : SampleList{iLevel}.nSource
    TempS{iSource, iLevel} = zeros(size(TempT{iSource, iLevel}));
  end

  for iSource = 1 : SampleList{iLevel}.nSource
    for iTarget = 1 : SampleList{iLevel}.nTarget
      iRKBlock = SampleList{iLevel}.RKIndex(iTarget, iSource );
      if( SampleList{iLevel}.SymFlag(iTarget, iSource ) == DIRECT )
	TargetInd = RKTree{iLevel}.RKTargetInd(iRKBlock);
	TempS{TargetInd, iLevel} = TempS{TargetInd, iLevel} + ...
	  RKTree{iLevel}.RKBlockH2M{iRKBlock} * TempT{iSource, iLevel};
      else
	TargetInd = RKTree{iLevel}.RKSourceInd(iRKBlock);
	TempS{TargetInd, iLevel} = TempS{TargetInd, iLevel} + ...
	  transpose(RKTree{iLevel}.RKBlockH2M{iRKBlock}) * ...
	  TempT{iSource, iLevel};
      end
    end %iTarget
  end % iSource

end

for iLevel = 1 : nLevelCur -2 
  for iSource = 1:SampleList{iLevel}.nSource
    for iChild = 1:SampleList{iLevel}.nChild
      TempS{SampleList{iLevel}.Child(iChild, iSource), iLevel+1} = ...
          TempS{SampleList{iLevel}.Child(iChild, iSource), iLevel+1} ...
          + RKTree{iLevel}.RKBlockH2T{iChild, iSource} * TempS{iSource, ...
                    iLevel};
    
    end
  end 
end %iLevel

iLevel = nLevelCur - 1;
nWidth = RKTree{iLevel}.nWidth;

% iRow is used because here it does not represent the column (source) direction
% of the matrix, 
for iRow = 1 : SampleList{iLevel}.nSource
  Row = SampleList{iLevel}.Source( :, iRow );
  TempR = RKTree{iLevel}.RKBlockH2U{iRow} * TempS{iRow,iLevel};
  v( 1+(Row(1)-1)*nWidth : Row(1)*nWidth, ...
    1+(Row(2)-1)*nWidth : Row(2)*nWidth, : ) = ...
    reshape( TempR, nWidth, nWidth, nSample );
end % iRow


