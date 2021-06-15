% TESTGREEN is a test subroutine.  It performs a block-by-block
% comparison between the Green's function of operator
%   H = -Delta + 1
% and the Green's function after uniform H^1 or H^2 compression.
%
% Adjustable Parameter:
%
% HMode              :    'H2'-> H^2 compression.
%                         'UNIFORM'-> uniform H^1 compression.
% N                  :    number of grid points per dimension.
% Nlevel             :    number of low rank level.


% Parameter setup
HMode = 'H2';

fprintf(' HMode = %s\n', HMode);

% Main test code 

if(1)
  N = 32;

  h = (1/N);
  V = 1 + rand(N,N);
  % V = ones(N,N);

  global Tree;
  global mvcount
  disp('Factorization');
  tic
  Tree = five_setup(V,N,h);
  toc
  disp('Extraction');
  tic
  DiagExact = five_extract(N,Tree);
  DiagExact = reshape( DiagExact, N*N, 1 );
  toc
end

if(1)
  disp('Generating Hamiltonian matrix');
  tic
  % Dirichlet boundary condition
  H = 1/h^2 * delsq(numgrid('S',N+2));
  % Add Periodic boundary condition
  j = 1;
  for i = 1 : N
    pos = N*(j-1) + i;
    pos1 = N * (N-1) + i;
    H(pos, pos1) = -1/h^2;
    H(pos1, pos) = -1/h^2;
  end
  i = 1;
  for j = 1 : N
    pos = N*(j-1) + i;
    pos1 = N * (j-1) + N;
    H(pos, pos1) = -1/h^2;
    H(pos1, pos) = -1/h^2;
  end
  
  % Add Potential
  H = H + spdiags(reshape(V,[],1), 0, N*N, N*N);
  toc
  
  % Inverse
  disp('Calculating inverse');
  tic
  G = inv(H);
  toc
end

if(1)
  svdfact = 1;
  epsilon = 1e-6;
  nLevel = 3;


  disp('Constructing H matrix');
  [DiagVec, DiagBlock, nDiagBlock, SampleList, RKTree] = ...
    ConstructHMatrix( N, nLevel, epsilon, svdfact );
end

% Compare the L^2 error for each block on each level.
if(1)
  MeshInd = reshape(1:N*N,N,N);
  GabserrL2 = cell(nLevel+1,1);
  GrelerrL2 = cell(nLevel+1,1);
  GabserrMax = cell(nLevel+1,1);
  GrelerrMax = cell(nLevel+1,1);
  GabserrL1 = cell(nLevel+1,1);
  GrelerrL1 = cell(nLevel+1,1);
  for iLevel = 1 : nLevel+1
    if( iLevel <= nLevel )
      GabserrL2{iLevel} = zeros(RKTree{iLevel}.nRKBlock,1);
      GrelerrL2{iLevel} = zeros(RKTree{iLevel}.nRKBlock,1);
      GabserrMax{iLevel} = zeros(RKTree{iLevel}.nRKBlock,1);
      GrelerrMax{iLevel} = zeros(RKTree{iLevel}.nRKBlock,1);
      GabserrL1{iLevel} = zeros(RKTree{iLevel}.nRKBlock,1);
      GrelerrL1{iLevel} = zeros(RKTree{iLevel}.nRKBlock,1);
      fprintf('Level %5i\n', iLevel);
      switch( HMode )
	case 'UNIFORM'
	  for iRKBlock = 1: RKTree{iLevel}.nRKBlock
	    Source = RKTree{iLevel}.RKSource(:,:,iRKBlock);  
	    Target = RKTree{iLevel}.RKTarget(:,:,iRKBlock);  
	    SourcePos = ...
	      reshape(MeshInd(Source(1,1):Source(1,2), ...
	      Source(2,1):Source(2,2)), [], 1);
	    TargetPos = ...
	      reshape(MeshInd(Target(1,1):Target(1,2), ...
	      Target(2,1):Target(2,2)), [], 1);
	    SourceInd = RKTree{iLevel}.RKSourceInd(iRKBlock);
	    TargetInd = RKTree{iLevel}.RKTargetInd(iRKBlock);
	    Gexact = G(TargetPos, SourcePos);
	    Gappro = RKTree{iLevel}.RKBlockUniU{TargetInd} * ...
	      RKTree{iLevel}.RKBlockUniM{iRKBlock} * ...
	      transpose(RKTree{iLevel}.RKBlockUniU{SourceInd});
	    GabserrL2{iLevel}(iRKBlock) = norm(Gexact - Gappro, 2);
	    GrelerrL2{iLevel}(iRKBlock) = GabserrL2{iLevel}(iRKBlock) ./ ...
	      norm(full(Gexact),2);
	    GabserrMax{iLevel}(iRKBlock) = max(max(abs(Gexact - Gappro)));
	    GrelerrMax{iLevel}(iRKBlock) = GabserrMax{iLevel}(iRKBlock) ./ ...
	      max(max(abs(Gexact)));
	    GabserrL1{iLevel}(iRKBlock) = sum(sum(abs(Gexact - Gappro))) ./ ...
	      numel(Gexact);
	    GrelerrL1{iLevel}(iRKBlock) = sum(sum(abs(Gexact - Gappro))) ./ ...
	      sum(sum(abs(Gexact)));
	    % VecGerr = sort(reshape(abs(Gexact - Gappro),[],1),'descend') ./ ...
	    % max(reshape(abs(Gexact-Gappro),[],1));
	    % VecGerr(100)
	    % pause
	  end % iRKBlock
	case 'H2'
	  for iRKBlock = 1: RKTree{iLevel}.nRKBlock
	    Source = RKTree{iLevel}.RKSource(:,:,iRKBlock);  
	    Target = RKTree{iLevel}.RKTarget(:,:,iRKBlock);  
	    SourcePos = ...
	      reshape(MeshInd(Source(1,1):Source(1,2), ...
	      Source(2,1):Source(2,2)), [], 1);
	    TargetPos = ...
	      reshape(MeshInd(Target(1,1):Target(1,2), ...
	      Target(2,1):Target(2,2)), [], 1);
	    SourceInd = RKTree{iLevel}.RKSourceInd(iRKBlock);
	    TargetInd = RKTree{iLevel}.RKTargetInd(iRKBlock);
	    Gexact = G(TargetPos, SourcePos);
	    Gappro = RKTree{iLevel}.RKBlockH2U{TargetInd} * ...
	      RKTree{iLevel}.RKBlockH2M{iRKBlock} * ...
	      transpose(RKTree{iLevel}.RKBlockH2U{SourceInd});
	    GabserrL2{iLevel}(iRKBlock) = norm(Gexact - Gappro, 2);
	    GrelerrL2{iLevel}(iRKBlock) = GabserrL2{iLevel}(iRKBlock) ./ ...
	      norm(full(Gexact),2);
	    GabserrMax{iLevel}(iRKBlock) = max(max(abs(Gexact - Gappro)));
	    GrelerrMax{iLevel}(iRKBlock) = GabserrMax{iLevel}(iRKBlock) ./ ...
	      max(max(abs(Gexact)));
	    GabserrL1{iLevel}(iRKBlock) = sum(sum(abs(Gexact - Gappro))) ./ ...
	      numel(Gexact);
	    GrelerrL1{iLevel}(iRKBlock) = sum(sum(abs(Gexact - Gappro))) ./ ...
	      sum(sum(abs(Gexact)));
	    % VecGerr = sort(reshape(abs(Gexact - Gappro),[],1),'descend') ./ ...
	    % max(reshape(abs(Gexact-Gappro),[],1));
	    % VecGerr(100)
	    % pause
	  end % iRKBlock
      end
    else
      GabserrL2{iLevel} = zeros(nDiagBlock,1);
      GrelerrL2{iLevel} = zeros(nDiagBlock,1);
      GabserrMax{iLevel} = zeros(nDiagBlock,1);
      GrelerrMax{iLevel} = zeros(nDiagBlock,1);
      GabserrL1{iLevel} = zeros(nDiagBlock,1);
      GrelerrL1{iLevel} = zeros(nDiagBlock,1);
      fprintf('DiagBlock\n');
      for iBlock = 1: nDiagBlock
	SourcePos = DiagBlock{iBlock}.indSource;
	TargetPos = DiagBlock{iBlock}.indTarget;
	Gexact = G(TargetPos, SourcePos);
	Gappro = DiagBlock{iBlock}.Mat;
	GabserrL2{iLevel}(iBlock) = norm(Gexact - Gappro, 2);
	GrelerrL2{iLevel}(iBlock) = GabserrL2{iLevel}(iBlock) ./ ...
	  norm(full(Gexact),2);
	GabserrMax{iLevel}(iBlock) = max(max(abs(Gexact - Gappro)));
	GrelerrMax{iLevel}(iBlock) = GabserrMax{iLevel}(iBlock) ./ ...
	  max(max(abs(Gexact)));
	GabserrL1{iLevel}(iBlock) = sum(sum(abs(Gexact - Gappro))) ./ ...
	  numel(Gexact);
	GrelerrL1{iLevel}(iBlock) = sum(sum(abs(Gexact - Gappro))) ./ ...
	  sum(sum(abs(Gexact)));
	% VecGerr = sort(reshape(abs(Gexact - Gappro),[],1),'descend') ./ ...
	% max(reshape(abs(Gexact-Gappro),[],1));
	% VecGerr(100)
	% pause
      end % iRKBlock
    end
  end % iLevel
    for iLevel = 1 : nLevel+1
      if( iLevel <= nLevel)
	disp(['Max absolute L^2 error of blocks at level ' num2str(iLevel)])
	max(GabserrL2{iLevel})
	disp(['Max relative L^2 error of blocks at level ' num2str(iLevel)])
	max(GrelerrL2{iLevel})
	disp(['Max absolute Max error of blocks at level ' num2str(iLevel)])
	max(GabserrMax{iLevel})
	disp(['Max relative Max error of blocks at level ' num2str(iLevel)])
	max(GrelerrMax{iLevel})
	disp(['Max absolute L^1 error of blocks at level ' num2str(iLevel)])
	max(GabserrL1{iLevel})
	disp(['Max relative L^1 error of blocks at level ' num2str(iLevel)])
	max(GrelerrL1{iLevel})
    else
      disp('Max absolute L^2 error of blocks at DiagBlock')
      max(GabserrL2{iLevel})
      disp('Max relative L^2 error of blocks at DiagBlock')
      max(GrelerrL2{iLevel})
      disp('Max absolute Max error of blocks at DiagBlock')
      max(GabserrMax{iLevel})
      disp('Max relative Max error of blocks at DiagBlock')
      max(GrelerrMax{iLevel})
      disp('Max absolute L^1 error of blocks at DiagBlock')
      max(GabserrL1{iLevel})
      disp('Max relative L^1 error of blocks at DiagBlock')
      max(GrelerrL1{iLevel})
    end

  end
  disp('Type "RKTree{i}.RKRankUniU" to see the rank distribution of i-th level after uniformization')
  fprintf('\n');
end

% Construct the low rank Green's function using brute force
if(0)
  fprintf('Constructing the H-matrix form Green''s function in a brute force way and calculate the 2-norm error \n');
  nSample = N*N;
  RMult = reshape(eye(N*N), N,N,N*N);
  Gbrute = WholeTreeApply( ...
    SampleList, RKTree, DiagBlock, RMult, N, nSample, nLevel+1, nDiagBlock);
  disp('absolute L^2 error')
  norm(Gbrute-G,2)
  disp('relative L^2 error')
  norm(Gbrute-G,2)./norm(full(G),2)
end
