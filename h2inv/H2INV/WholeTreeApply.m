function v = WholeTreeApply( ...
  SampleList, RKTree, DiagBlock, u, N, nSample, nLevelCur, nDiagBlock )
% WholeTreeApply is the complete Matrix vector multiplication, including
% both the contribution of low rank matrices and the contribution from
% the diagonal and  nearest off-diagonal matrices. 
% ********************************************************
% INPUT PARAMETERS
% ********************************************************
% SampleList:   Cell structure, containing the information of the
%               decomposd domain.
% RKTree    :   Cell structure, containing the low rank matrices and
%               related information. 
% DiagBlock :   Cell structure, containing the part of Green's function
%               that is not represented by low rank matrics. This part
%               includes the diagonal blocks as well as off-diagonal
%               blocks.
% u         :   N * N * nSample array (3D). Vectors to be multiplied.
%               Interpreted as sources.
% nSample   :   Number of right hand sides.
% nLevelCur :   The length of RKTree + 1.
% nDiagBlock:   The length of DiagBlock.
% ********************************************************
% OUTPUT PARAMETERS
% ********************************************************
% v         :   (N * N) * nSample array (2D). Vectors after
%               matrix-vector multiplcation.  Interpreted as targets.

v = RKTreeApplyH2( SampleList, RKTree, u, N, nLevelCur, nSample );
v = reshape(v, N*N, nSample);
w = reshape(u, N*N, nSample);
for i = 1 : nDiagBlock
  v(DiagBlock{i}.indTarget,:) = v(DiagBlock{i}.indTarget,:) + ...
    DiagBlock{i}.Mat * w(DiagBlock{i}.indSource,:);
end
