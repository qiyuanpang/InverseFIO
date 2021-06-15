% TESTTRANSFER is a test subroutine.  By running testGreen.m with 
% nLevel = 3, The matrices
%   H2U * H2M * transpose(H2U) at level 2
% is compared with
%   H2U * H2T * H2M * transpose(H2U * H2T) at level 3.
% Also
%   H2U * H2M * transpose(H2U) at level 1
% is compared with
%   H2U * H2T * H2T * H2M * transpose(H2U * H2T * H2T) at level 3.
%
% The comparison is done in a block-by-block fashion.  If there is no
% output, it means the transfer matrices H2T are calculated properly.


%% Level 2 test

if(0)
  disp(' Level 2 test' );
  for iRKBlock = 1 : RKTree{2}.nRKBlock
    SourceInd = RKTree{2}.RKSourceInd(iRKBlock);
    TargetInd = RKTree{2}.RKTargetInd(iRKBlock);

    TA1 = cell(SampleList{2}.nChild);
    TA2 = cell(SampleList{2}.nChild);

    nWidth      = RKTree{2}.nWidth;
    nWidthChild = RKTree{3}.nWidth;
    tmpG = reshape((1:nWidth*nWidth), nWidth, nWidth);
    IND = cell(SampleList{2}.nChild);
    IND{1} = reshape(tmpG(1:nWidthChild, 1:nWidthChild),[],1);
    IND{2} = reshape(tmpG(nWidthChild+1:2*nWidthChild, 1:nWidthChild),[],1);
    IND{3} = reshape(tmpG(1:nWidthChild, nWidthChild+1:2*nWidthChild),[],1);
    IND{4} = reshape(tmpG(nWidthChild+1:2*nWidthChild, nWidthChild+1:2*nWidthChild),[],1);

    for iChild = 1 : SampleList{2}.nChild
      ChildTargetInd = SampleList{2}.Child(iChild, TargetInd);
      TA1{iChild} = RKTree{3}.RKBlockH2U{ChildTargetInd} * ...
	transpose(RKTree{3}.RKBlockH2U{ChildTargetInd}) * ...
	RKTree{2}.RKBlockH2UChild{iChild,TargetInd};
    end
    for iChild = 1 : SampleList{2}.nChild
      ChildSourceInd = SampleList{2}.Child(iChild, SourceInd);
      TA2{iChild} = RKTree{3}.RKBlockH2U{ChildSourceInd} * ...
	transpose(RKTree{3}.RKBlockH2U{ChildSourceInd}) * ...
	RKTree{2}.RKBlockH2UChild{iChild,SourceInd};
    end

    % TB1 = RKTree{3}.RKBlockH2U{1} * RKTree{2}.RKBlockH2T{1,1};
    % TB2 = RKTree{3}.RKBlockH2U{5} * RKTree{2}.RKBlockH2T{1,3};
    % 
    % B = TB1 * RKTree{2}.RKBlockH2M{indM} * transpose(TB2)

    C = RKTree{2}.RKBlockH2U{TargetInd} * RKTree{2}.RKBlockH2M{iRKBlock} * ...
      transpose(RKTree{2}.RKBlockH2U{SourceInd});



    for iChild = 1 : SampleList{2}.nChild
      for jChild = 1 : SampleList{2}.nChild
	A = TA1{iChild} * RKTree{2}.RKBlockH2M{iRKBlock} * ...
	  transpose(TA2{jChild});
	E = C(IND{iChild},IND{jChild})-A;
	if( norm(E) > 1e-10 )
	  iChild,jChild, E
	end
      end
    end

  end
end

%% Level 1 test

if(1) 
  disp(' Level 1 test' );
  for iRKBlock = 1 : RKTree{2}.nRKBlock
    SourceInd = RKTree{2}.RKSourceInd(iRKBlock);
    TargetInd = RKTree{2}.RKTargetInd(iRKBlock);

    TA1 = cell(SampleList{2}.nChild);
    TA2 = cell(SampleList{2}.nChild);

    nWidth      = RKTree{2}.nWidth;
    nWidthChild = RKTree{3}.nWidth;
    tmpG = reshape((1:nWidth*nWidth), nWidth, nWidth);
    IND = cell(SampleList{2}.nChild);
    IND{1} = reshape(tmpG(1:nWidthChild, 1:nWidthChild),[],1);
    IND{2} = reshape(tmpG(nWidthChild+1:2*nWidthChild, 1:nWidthChild),[],1);
    IND{3} = reshape(tmpG(1:nWidthChild, nWidthChild+1:2*nWidthChild),[],1);
    IND{4} = reshape(tmpG(nWidthChild+1:2*nWidthChild, nWidthChild+1:2*nWidthChild),[],1);




    for iChild = 1 : SampleList{2}.nChild
      ChildTargetInd = SampleList{2}.Child(iChild, TargetInd);
      TA1{iChild} = RKTree{3}.RKBlockH2U{ChildTargetInd} * ...
	transpose(RKTree{3}.RKBlockH2U{ChildTargetInd}) * ...
	RKTree{2}.RKBlockH2UChild{iChild,TargetInd};
    end
    for iChild = 1 : SampleList{2}.nChild
      ChildSourceInd = SampleList{2}.Child(iChild, SourceInd);
      TA2{iChild} = RKTree{3}.RKBlockH2U{ChildSourceInd} * ...
	transpose(RKTree{3}.RKBlockH2U{ChildSourceInd}) * ...
	RKTree{2}.RKBlockH2UChild{iChild,SourceInd};
    end

    % TB1 = RKTree{3}.RKBlockH2U{1} * RKTree{2}.RKBlockH2T{1,1};
    % TB2 = RKTree{3}.RKBlockH2U{5} * RKTree{2}.RKBlockH2T{1,3};
    % 
    % B = TB1 * RKTree{2}.RKBlockH2M{indM} * transpose(TB2)

    C = RKTree{2}.RKBlockH2U{TargetInd} * RKTree{2}.RKBlockH2M{iRKBlock} * ...
      transpose(RKTree{2}.RKBlockH2U{SourceInd});



    for iChild = 1 : SampleList{2}.nChild
      for jChild = 1 : SampleList{2}.nChild
	A = TA1{iChild} * RKTree{2}.RKBlockH2M{iRKBlock} * ...
	  transpose(TA2{jChild});
	E = C(IND{iChild},IND{jChild})-A;
	if( norm(E) > 1e-10 )
	  iChild,jChild, E
	end
      end
    end

  end
end

%% Level 1 test

if(1) 
  for iRKBlock = 1 : RKTree{1}.nRKBlock
    SourceInd = RKTree{1}.RKSourceInd(iRKBlock);
    TargetInd = RKTree{1}.RKTargetInd(iRKBlock);

    nWidth1 = RKTree{1}.nWidth;
    nWidth2 = RKTree{2}.nWidth;
    nWidth3 = RKTree{3}.nWidth;

    tmpG = reshape((1:nWidth1*nWidth1), nWidth1, nWidth1);
    IND1 = cell(4,1);
    IND1{1} = reshape(tmpG(1:nWidth2, 1:nWidth2),[],1);
    IND1{2} = reshape(tmpG(nWidth2+1:2*nWidth2, 1:nWidth2),[],1);
    IND1{3} = reshape(tmpG(1:nWidth2, nWidth2+1:2*nWidth2),[],1);
    IND1{4} = reshape(tmpG(nWidth2+1:2*nWidth2, nWidth2+1:2*nWidth2),[],1);

    tmpG = reshape((1:nWidth2*nWidth2), nWidth2, nWidth2);
    IND2 = cell(4,1);
    IND2{1} = reshape(tmpG(1:nWidth3, 1:nWidth3),[],1);
    IND2{2} = reshape(tmpG(nWidth3+1:2*nWidth3, 1:nWidth3),[],1);
    IND2{3} = reshape(tmpG(1:nWidth3, nWidth3+1:2*nWidth3),[],1);
    IND2{4} = reshape(tmpG(nWidth3+1:2*nWidth3, nWidth3+1:2*nWidth3),[],1);




    for iChild1 = 1 : 4
      for jChild1 = 1 : 4
	for iChild2 = 1 : 4
	  for jChild2 = 1 : 4
	    ChildTargetInd1 = SampleList{1}.Child(iChild1, TargetInd);
	    ChildTargetInd2 = SampleList{2}.Child(iChild2, ChildTargetInd1);
	    ChildSourceInd1 = SampleList{1}.Child(jChild1, SourceInd);
	    ChildSourceInd2 = SampleList{2}.Child(jChild2, ChildSourceInd1);

	    % TA1 = ...
	      % RKTree{3}.RKBlockH2U{ChildTargetInd2} * ...
	      % transpose(RKTree{3}.RKBlockH2U{ChildTargetInd2}) * ...
	      % RKTree{2}.RKBlockH2UChild{iChild2,ChildTargetInd1} * ...
	      % transpose(RKTree{2}.RKBlockH2U{ChildTargetInd1}) * ...
	      % RKTree{1}.RKBlockH2UChild{iChild1,TargetInd};
	    % TA2 = ...
	      % RKTree{3}.RKBlockH2U{ChildSourceInd2} * ...
	      % transpose(RKTree{3}.RKBlockH2U{ChildSourceInd2}) * ...
	      % RKTree{2}.RKBlockH2UChild{jChild2,ChildSourceInd1} * ...
	      % transpose(RKTree{2}.RKBlockH2U{ChildSourceInd1}) * ...
	      % RKTree{1}.RKBlockH2UChild{jChild1,SourceInd};
	    % A = TA1 * RKTree{1}.RKBlockH2M{iRKBlock} * ...
	      % transpose(TA2);

	    TB1 = RKTree{3}.RKBlockH2U{ChildTargetInd2} * ...
	      RKTree{2}.RKBlockH2T{iChild2, ChildTargetInd1} * ...
	      RKTree{1}.RKBlockH2T{iChild1, TargetInd}; 
	    TB2 = RKTree{3}.RKBlockH2U{ChildSourceInd2} * ...
	      RKTree{2}.RKBlockH2T{jChild2, ChildSourceInd1} * ...
	      RKTree{1}.RKBlockH2T{jChild1, SourceInd}; 
	    B = TB1 * RKTree{1}.RKBlockH2M{iRKBlock} * ...
	      transpose(TB2);

	    C = RKTree{1}.RKBlockH2U{TargetInd} * RKTree{1}.RKBlockH2M{iRKBlock} * ...
	      transpose(RKTree{1}.RKBlockH2U{SourceInd});
	    Tmp = C(IND1{iChild1}, IND1{jChild1});

	    E = Tmp(IND2{iChild2},IND2{jChild2})-B;
	    % if( norm(E) > 1e-10 )
	      iChild1,jChild1, iChild2, jChild2,E, Tmp(IND2{iChild2},IND2{jChild2})
	      pause
	    % end
	  end
	end
      end
    end


  end
end
