% Covariance matrix on the unit line, squared exponential kernel.
%
% This example considers the covariance matrix between equispaced points on the
% unit line with squared exponential covariance function. The matrix is square,
% real, positive definite, and Toeplitz.
%
% We remark on a few key departures from the PDE setting:
%
%   - the covariance kernel has a length scale parameter
%   - the covariance matrix has a "nugget" (identity perturbation) to improve
%       conditioning
%   - the proxy points sample a local corona since there is no Green's theorem
%       (but the kernel is real-analytic)
%   - computing explicit entries of the inverse covariance matrix (i.e.,
%       precision matrix) is of interest
%
% This demo does the following in order:
%
%   - compress the matrix
%   - check multiply error/time
%   - build/factor extended sparsification
%   - check solve error/time using extended sparsification
%   - do selected inversion
%   - do diagonal inversion (i.e., selected inversion for the diagonal)
%
% Inputs (defaults are used if not provided or set empty):
%
%   - N: number of points (default: N = 16384)
%   - OCC: tree occupancy parameter (default: OCC = 128)
%   - P: half-number of proxy points (default: P = 8)
%   - RANK_OR_TOL: local precision parameter (default: RANK_OR_TOL = 1e-12)
%   - TMAX: ID interpolation matrix entry bound (default: TMAX = 2)
%   - SYMM: symmetry parameter (default: SYMM = 'P')
%   - NOISE: nugget effect (default: NOISE = 1e-2)
%   - SCALE: kernel length scale (default: SCALE = 100)

function cov_line1(N,occ,p,rank_or_tol,Tmax,symm,noise,scale)

  % set default parameters
  if nargin < 1 || isempty(N), N = 16384; end
  if nargin < 2 || isempty(occ), occ = 128; end
  if nargin < 3 || isempty(p), p = 8; end
  if nargin < 4 || isempty(rank_or_tol), rank_or_tol = 1e-12; end
  if nargin < 5 || isempty(Tmax), Tmax = 2; end
  if nargin < 6 || isempty(symm), symm = 'p'; end
  if nargin < 7 || isempty(noise), noise = 1e-2; end
  if nargin < 8 || isempty(scale), scale = 100; end

  % initialize
  x = (1:N)/N;                                          % grid points
  proxy = linspace(1.5,2.5,p); proxy = [-proxy proxy];  % proxy points
  % reference proxy points are for unit box [-1, 1]

  % compress matrix
  Afun = @(i,j)Afun_(i,j,x,noise,scale);
  pxyfun = @(rc,rx,cx,slf,nbr,l,ctr)pxyfun_(rc,rx,cx,slf,nbr,l,ctr,proxy,scale);
  opts = struct('Tmax',Tmax,'symm',symm,'verb',1);
  tic; F = rskel(Afun,x,x,occ,rank_or_tol,pxyfun,opts); t = toc;
  w = whos('F'); mem = w.bytes/1e6;
  fprintf('rskel time/mem: %10.4e (s) / %6.2f (MB)\n',t,mem)

  % set up reference FFT multiplication
  a = Afun(1:N,1);
  B = zeros(2*N-1,1);  % zero-pad
  B(1:N) = a;
  B(N+1:end) = flipud(a(2:N));
  G = fft(B);
  mv = @(x)mv_(G,x);

  % test accuracy using randomized power method
  X = rand(N,1);
  X = X/norm(X);

  % NORM(A - F)/NORM(A)
  tic; rskel_mv(F,X); t = toc;  % for timing
  err = snorm(N,@(x)(mv(x) - rskel_mv(F,x)),[],[],1);
  err = err/snorm(N,mv,[],[],1);
  fprintf('rskel_mv err/time: %10.4e / %10.4e (s)\n',err,t)

  % build extended sparsification
  tic; [A,p,q] = rskel_xsp(F); t = toc;
  w = whos('A'); mem = w.bytes/1e6;
  fprintf('rskel_xsp:\n')
  fprintf('  build time/mem: %10.4e (s) / %6.2f (MB)\n',t,mem);

  % factor extended sparsification
  dolu = F.symm == 'n';  % LU or LDL?
  if ~dolu && isoctave()
    warning('No LDL in Octave; using LU.')
    dolu = 1;
    A = A + tril(A,-1)';
  end
  FA = struct('p',p,'q',q,'lu',dolu);
  tic
  if dolu, [FA.L,FA.U,FA.P] = lu(A);
  else,    [FA.L,FA.D,FA.P] = ldl(A);
  end
  t = toc;
  w = whos('FA'); mem = w.bytes/1e6;
  fprintf('  factor time/mem: %10.4e (s) / %6.2f (MB)\n',t,mem)
  sv = @(x,trans)sv_(FA,x,trans);  % linear solve function

  % NORM(INV(A) - INV(F))/NORM(INV(A)) <= NORM(I - A*INV(F))
  tic; sv(X,'n'); t = toc;  % for timing
  err = snorm(N,@(x)(x - mv(sv(x,'n'))),@(x)(x - sv(mv(x),'c')));
  fprintf('  solve err/time: %10.4e / %10.4e (s)\n',err,t)

  % prepare for selected/diagonal inversion
  % note: reported error is with respect to compressed solve -- should ideally
  %       be close O(eps), independently of compression tolerance
  m = 128;           % number of entries
  r = randi(N,m,2);  % random row/col index for each entry
  % reference comparison from compressed solve against coordinate vectors
  X = zeros(N,m);
  for i = 1:m, X(r(i,2),i) = 1; end
  Y = sv(X,'n');  % subselect r(:,1) or r(:,2) rows for selected/diagonal inv
  % coordinate vector storage in sparse form
  ei = sparse([],[],[],size(A,1),1,1);
  ej = sparse([],[],[],size(A,1),1,1);
  % solution storage
  S = zeros(m,1);
  T = zeros(m,1);

  % selected inversion
  % algorithm: sparse dot product of (inverse) matrix factor row/cols
  % complexity: O(LOG(N)) for each entry
  % note: error can sometimes be poor since exact entries are small
  tic
  for i = 1:m
    ei(FA.q(r(i,1))) = 1;
    ej(FA.p(r(i,2))) = 1;
    if dolu, S(i) = (FA.U'\ei)'*(FA.L\(FA.P*ej));
    else,    S(i) = (FA.L\(FA.P\ei))'*(FA.D\(FA.L\(FA.P\ej)));
    end
    ei(FA.q(r(i,1))) = 0;
    ej(FA.p(r(i,2))) = 0;
  end
  t = toc/m;  % average time
  for i = 1:m, T(i) = Y(r(i,1),i); end
  err = norm(S - T)/norm(T);
  fprintf('selinv err/avg time: %10.4e / %10.4e (s)\n',err,t)

  % diagonal inversion -- same as selinv but for most "important" entries
  tic
  for i = 1:m
    ei(FA.p(r(i,2))) = 1;
    if dolu, S(i) = (FA.U'\ei)'*(FA.L\(FA.P*ei));
    else,    S(i) = (FA.L\(FA.P\ei))'*(FA.D\(FA.L\(FA.P\ei)));
    end
    ei(FA.p(r(i,2))) = 0;
  end
  t = toc/m;  % average time
  for i = 1:m, T(i) = Y(r(i,2),i); end
  err = norm(S - T)/norm(T);
  fprintf('diaginv err/avg time: %10.4e / %10.4e (s)\n',err,t)
end

% kernel function
function K = Kfun(x,y,scale)
  dr = scale*abs(x' - y);  % scaled distance
  K = exp(-0.5*dr.^2);
end

% matrix entries
function A = Afun_(i,j,x,noise,scale)
  A = Kfun(x(:,i),x(:,j),scale);
  [I,J] = ndgrid(i,j);
  idx = I == J;
  A(idx) = A(idx) + noise^2;  % modify diagonal with "nugget"
end

% proxy function
function [Kpxy,nbr] = pxyfun_(rc,rx,cx,slf,nbr,l,ctr,proxy,scale)
  pxy = proxy.*l + ctr;  % scale and translate reference points
  N = size(rx,2);
  if rc == 'r'
    Kpxy = Kfun(rx(:,slf),pxy,scale);
    dr = cx(nbr) - ctr;
  else
    Kpxy = Kfun(pxy,cx(:,slf),scale);
    dr = rx(nbr) - ctr;
  end
  % proxy points form interval of scaled radius 1.5 around current box
  % keep among neighbors only those within interval
  nbr = nbr(abs(dr)/l < 1.5);
end

% FFT multiplication
function y = mv_(F,x)
  n = length(x);
  y = ifft(F.*fft(x,2*n-1));
  y = y(1:n);
end

% sparse LU/LDL solve
function Y = sv_(F,X,trans)
  N = size(X,1);
  if trans == 'n', p = F.p; q = F.q;
  else,            p = F.q; q = F.p;
  end
  X = [X(p,:); zeros(size(F.L,1)-N,size(X,2))];
  if F.lu
    if trans == 'n', Y = F.U \(F.L \(F.P *X));
    else,            Y = F.P'*(F.L'\(F.U'\X));
    end
  else
    Y = F.P*(F.L'\(F.D\(F.L\(F.P'*X))));
  end
  Y = Y(1:N,:);
  Y(q,:) = Y;
end