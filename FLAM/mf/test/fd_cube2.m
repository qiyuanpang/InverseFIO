% Seven-point stencil on the unit cube, variable-coefficient Poisson equation,
% Dirichlet boundary conditions.
%
% This is basically the 3D analogue of FD_SQUARE2.
%
% Inputs (defaults are used if not provided or set empty):
%
%   - N: number of points + 1 in each dimension (default: N = 32)
%   - OCC: tree occupancy parameter (default: OCC = 4)
%   - SYMM: symmetry parameter (default: SYMM = 'P')
%   - DOITER: whether to run unpreconditioned CG (default: DOITER = 1)
%   - DIAGMODE: diagonal extraction mode - 0: skip; 1: matrix unfolding; 2:
%       sparse apply/solves (default: DIAGMODE = 0)

function fd_cube2(n,occ,symm,doiter,diagmode)

  % set default parameters
  if nargin < 1 || isempty(n), n = 32; end
  if nargin < 2 || isempty(occ), occ = 4; end
  if nargin < 3 || isempty(symm), symm = 'p'; end
  if nargin < 4 || isempty(doiter), doiter = 1; end
  if nargin < 5 || isempty(diagmode), diagmode = 0; end

  % initialize
  N = (n - 1)^3;  % total number of grid points

  % set up conductivity field
  a = zeros(n+1,n+1,n+1);
  A = rand(n-1,n-1,n-1);  % random field
  A = fftn(A,[2*n-3 2*n-3 2*n-3]);
  [X,Y,Z] = ndgrid(0:n-2);
  % Gaussian smoothing over 4 grid points
  C = gausspdf(X,0,4).*gausspdf(Y,0,4).*gausspdf(Z,0,4);
  B = zeros(2*n-3,2*n-3,2*n-3);
  B(1:n-1,1:n-1,1:n-1) = C;
  B(1:n-1,1:n-1,n:end) = C( :   , :   ,2:n-1);
  B(1:n-1,n:end,1:n-1) = C( :   ,2:n-1, :   );
  B(1:n-1,n:end,n:end) = C( :   ,2:n-1,2:n-1);
  B(n:end,1:n-1,1:n-1) = C(2:n-1, :   , :   );
  B(n:end,1:n-1,n:end) = C(2:n-1, :   ,2:n-1);
  B(n:end,n:end,1:n-1) = C(2:n-1,2:n-1, :   );
  B(n:end,n:end,n:end) = C(2:n-1,2:n-1,2:n-1);
  B(:,:,n:end) = flipdim(B(:,:,n:end),3);
  B(:,n:end,:) = flipdim(B(:,n:end,:),2);
  B(n:end,:,:) = flipdim(B(n:end,:,:),1);
  B = fftn(B);
  A = ifftn(A.*B);        % convolution in Fourier domain
  A = A(1:n-1,1:n-1,1:n-1);
  idx = A > median(A(:));
  A( idx) = 1e+2;         % set upper 50% to something large
  A(~idx) = 1e-2;         % set lower 50% to something small
  a(2:n,2:n,2:n) = A;
  clear X Y Z A B C

  % set up sparse matrix
  idx = zeros(n+1,n+1,n+1);  % index mapping to each point, including "ghosts"
  idx(2:n,2:n,2:n) = reshape(1:N,n-1,n-1,n-1);
  mid = 2:n;    % "middle" indices -- interaction with self
  lft = 1:n-1;  % "left"   indices -- interaction with one below
  rgt = 3:n+1;  % "right"  indices -- interaction with one above
  I = idx(mid,mid,mid);
  % interactions with ...
  Jl = idx(lft,mid,mid); Sl = -0.5*(a(lft,mid,mid) + a(mid,mid,mid));  % left
  Jr = idx(rgt,mid,mid); Sr = -0.5*(a(rgt,mid,mid) + a(mid,mid,mid));  % right
  Ju = idx(mid,lft,mid); Su = -0.5*(a(mid,lft,mid) + a(mid,mid,mid));  % up
  Jd = idx(mid,rgt,mid); Sd = -0.5*(a(mid,rgt,mid) + a(mid,mid,mid));  % down
  Jf = idx(mid,mid,lft); Sf = -0.5*(a(mid,mid,lft) + a(mid,mid,mid));  % front
  Jb = idx(mid,mid,rgt); Sb = -0.5*(a(mid,mid,rgt) + a(mid,mid,mid));  % back
  Jm = idx(mid,mid,mid); Sm = -(Sl + Sr + Sd + Su + Sb + Sf);  % middle (self)
  % combine all interactions
  I = [ I(:);  I(:);  I(:);  I(:);  I(:);  I(:);  I(:)];
  J = [Jl(:); Jr(:); Ju(:); Jd(:); Jf(:); Jb(:); Jm(:)];
  S = [Sl(:); Sr(:); Su(:); Sd(:); Sf(:); Sb(:); Sm(:)];
  % remove ghost interactions
  idx = find(J > 0); I = I(idx); J = J(idx); S = S(idx);
  A = sparse(I,J,S,N,N);
  clear idx Jl Sl Jr Sr Ju Su Jd Sd Jf Sf Jb Sb Im Sm I J S

  % factor matrix
  opts = struct('symm',symm,'verb',1);
  tic; F = mf3(A,n,occ,opts); t = toc;
  w = whos('F'); mem = w.bytes/1e6;
  fprintf('mf3 time/mem: %10.4e (s) / %6.2f (MB)\n',t,mem)

  % test accuracy using randomized power method
  X = rand(N,1);
  X = X/norm(X);

  % NORM(A - F)/NORM(A)
  tic; mf_mv(F,X); t = toc;  % for timing
  err = snorm(N,@(x)(A*x - mf_mv(F,x)),[],[],1);
  err = err/snorm(N,@(x)(A*x),[],[],1);
  fprintf('mf_mv: %10.4e / %10.4e (s)\n',err,t)

  % NORM(INV(A) - INV(F))/NORM(INV(A)) <= NORM(I - A*INV(F))
  tic; mf_sv(F,X); t = toc;  % for timing
  err = snorm(N,@(x)(x - A*mf_sv(F,x)),@(x)(x - mf_sv(F,A*x,'c')));
  fprintf('mf_sv: %10.4e / %10.4e (s)\n',err,t)

  % test Cholesky accuracy -- error is w.r.t. compressed apply/solve
  if strcmpi(symm,'p')
    % NORM(F - C*C')/NORM(F)
    tic; mf_cholmv(F,X); t = toc;  % for timing
    err = snorm(N,@(x)(mf_mv(F,x) - mf_cholmv(F,mf_cholmv(F,x,'c'))),[],[],1);
    err = err/snorm(N,@(x)mf_mv(F,x),[],[],1);
    fprintf('mf_cholmv: %10.4e / %10.4e (s)\n',err,t)

    % NORM(INV(F) - INV(C')*INV(C))/NORM(INV(F))
    tic; mf_cholsv(F,X); t = toc;  % for timing
    err = snorm(N,@(x)(mf_sv(F,x) - mf_cholsv(F,mf_cholsv(F,x),'c')),[],[],1);
    err = err/snorm(N,@(x)mf_sv(F,x),[],[],1);
    fprintf('mf_cholsv: %10.4e / %10.4e (s)\n',err,t)
  end

  % run unpreconditioned CG
  B = A*X;
  iter = nan;
  if doiter, [~,~,~,iter] = pcg(@(x)(A*x),B,1e-12,128); end

  % run preconditioned CG
  tic; [Y,~,~,piter] = pcg(@(x)(A*x),B,1e-12,32,@(x)mf_sv(F,x)); t = toc;
  err1 = norm(X - Y)/norm(X);
  err2 = norm(B - A*Y)/norm(B);
  fprintf('cg:\n')
  fprintf('  soln/resid err/time: %10.4e / %10.4e / %10.4e (s)\n', ...
          err1,err2,t)
  fprintf('  precon/unprecon iter: %d / %d\n',piter,iter)

  % compute log-determinant
  tic; ld = mf_logdet(F); t = toc;
  fprintf('mf_logdet:\n')
  fprintf('  real/imag: %22.16e / %22.16e\n',real(ld),imag(ld))
  fprintf('  time: %10.4e (s)\n',t)

  if diagmode > 0
    % prepare for diagonal extraction
    opts = struct('verb',1);
    m = min(N,128);  % number of entries to check against
    r = randperm(N); r = r(1:m);
    % reference comparison from compressed solve against coordinate vectors
    X = zeros(N,m);
    for i = 1:m, X(r(i),i) = 1; end
    E = zeros(m,1);  % solution storage
    if diagmode == 1, fprintf('mf_diag:\n')
    else,             fprintf('mf_spdiag:\n')
    end

    % extract diagonal
    tic;
    if diagmode == 1, D = mf_diag(F,0,opts);
    else,             D = mf_spdiag(F);
    end
    t = toc;
    Y = mf_mv(F,X);
    for i = 1:m, E(i) = Y(r(i),i); end
    err = norm(D(r) - E)/norm(E);
    fprintf('  fwd: %10.4e / %10.4e (s)\n',err,t)

    % extract diagonal of inverse
    tic;
    if diagmode == 1, D = mf_diag(F,1,opts);
    else,             D = mf_spdiag(F,1);
    end
    t = toc;
    Y = mf_sv(F,X);
    for i = 1:m, E(i) = Y(r(i),i); end
    err = norm(D(r) - E)/norm(E);
    fprintf('  inv: %10.4e / %10.4e (s)\n',err,t)
  end
end