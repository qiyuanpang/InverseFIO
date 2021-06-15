% GENKERNEL is a test subroutine.  After running testGreen.m one obtains
% SampleList and RKTree structure.  Then by changing iRKBlock and iLevel
% one can construct a Kernel that contains specifically one low-rank
% block for test purpose.

global Kernel;
Kernel = eye(N*N);
iRKBlock = 5;
iLevel   = 2;
iTargetArr = (RKTree{iLevel}.RKTarget(1,1,iRKBlock) : ...
  RKTree{iLevel}.RKTarget(1,2,iRKBlock))';
jTargetArr = RKTree{iLevel}.RKTarget(2,1,iRKBlock) : ...
  RKTree{iLevel}.RKTarget(2,2,iRKBlock);
indTarget = repmat( (jTargetArr-1)*N, length(iTargetArr), 1 ) + ...
  repmat( iTargetArr, 1, length(jTargetArr) );

iSourceArr = (RKTree{iLevel}.RKSource(1,1,iRKBlock) : ...
  RKTree{iLevel}.RKSource(1,2,iRKBlock))';
jSourceArr = RKTree{iLevel}.RKSource(2,1,iRKBlock) : ...
  RKTree{iLevel}.RKSource(2,2,iRKBlock);
indSource = repmat( (jSourceArr-1)*N, length(iSourceArr), 1 ) + ...
  repmat( iSourceArr, 1, length(jSourceArr) );


Kernel(indTarget(:), indSource(:)) = ones(numel(indTarget), ...
  numel(indSource));


Kernel = 0.5*(Kernel + Kernel');
