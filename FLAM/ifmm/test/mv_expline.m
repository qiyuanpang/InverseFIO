% Pairwise interactions on an exponentially graded line, Laplace kernel.
%
% This is basically the same as MV_LINE but now with points placed on a
% dyadically refined grid. The resulting geometry is highly non-uniform and
% provides a good test for adaptivity.
%
% Inputs (defaults are used if not provided or set empty):
%
%   - N: number of points (default: N = 64)
%   - OCC: tree occupancy parameter (default: OCC = 32)
%   - P: half-number of proxy points (default: P = 8)
%   - RANK_OR_TOL: local precision parameter (default: RANK_OR_TOL = 1e-12)
%   - TMAX: ID interpolation matrix entry bound (default: TMAX = 2)
%   - NEAR: near-field compression parameter (default: NEAR = 0)
%   - STORE: storage parameter (default: STORE = 'N')
%   - SYMM: symmetry parameter (default: SYMM = 'S')

function mv_expline(N,occ,p,rank_or_tol,Tmax,near,store,symm)

  % set default parameters
  if nargin < 1 || isempty(N), N = 64; end
  if nargin < 2 || isempty(occ), occ = 32; end
  if nargin < 3 || isempty(p), p = 8; end
  if nargin < 4 || isempty(rank_or_tol), rank_or_tol = 1e-12; end
  if nargin < 5 || isempty(Tmax), Tmax = 2; end
  if nargin < 6 || isempty(near), near = 0; end
  if nargin < 7 || isempty(store), store = 'n'; end
  if nargin < 8 || isempty(symm), symm = 's'; end

  % initialize
  x = 2.^(-(1:N));                                      % row/col points
  proxy = linspace(1.5,2.5,p); proxy = [-proxy proxy];  % proxy points
  % reference proxy points are for unit box [-1, 1]

  % compress matrix
  Afun = @(i,j)Afun_(i,j,x);
  pxyfun = @(rc,rx,cx,slf,nbr,l,ctr)pxyfun_(rc,rx,cx,slf,nbr,l,ctr,proxy);
  opts = struct('Tmax',Tmax,'near',near,'store',store,'symm',symm,'verb',1);
  tic; F = ifmm(Afun,x,x,occ,rank_or_tol,pxyfun,opts); t = toc;
  w = whos('F'); mem = w.bytes/1e6;
  fprintf('ifmm time/mem: %10.4e (s) / %6.2f (MB)\n',t,mem)

  % set up accuracy tests
  A = Afun(1:N,1:N);  % full matrix is okay since problem size is small
  X1 = rand(N,1); X1 = X1/norm(X1);  % for timing
  X = rand(N,16); X = X/norm(X);  % test against 16 vectors for robustness
  r = randperm(N); r = r(1:min(N,128));  % check up to 128 rows in result

  % test matrix apply accuracy
  tic; ifmm_mv(F,X1,Afun,'n'); t = toc;  % for timing
  Y = ifmm_mv(F,X,Afun,'n');
  Z = Afun(r,1:N)*X;
  err = norm(Z - Y(r,:))/norm(Z);
  fprintf('ifmm_mv:\n')
  fprintf('  multiply err/time: %10.4e / %10.4e (s)\n',err,t)

  % test matrix adjoint apply accuracy
  tic; ifmm_mv(F,X1,Afun,'c'); t = toc;  % for timing
  Y = ifmm_mv(F,X,Afun,'c');
  Z = Afun(1:N,r)'*X;
  err = norm(Z - Y(r,:))/norm(Z);
  fprintf('  adjoint multiply err/time: %10.4e / %10.4e (s)\n',err,t)
end

% kernel function
function K = Kfun(x,y)
  dr = abs(x' - y);
  K = dr;
end

% matrix entries
function A = Afun_(i,j,x)
  A = Kfun(x(:,i),x(:,j));
end

% proxy function
function [Kpxy,nbr] = pxyfun_(rc,rx,cx,slf,nbr,l,ctr,proxy)
  pxy = proxy.*l + ctr;  % scale and translate reference points
  if rc == 'r'
    Kpxy = Kfun(rx(:,slf),pxy);
    dr = cx(nbr) - ctr;
  else
    Kpxy = Kfun(pxy,cx(:,slf));
    dr = rx(nbr) - ctr;
  end
  % proxy points form interval of scaled radius 1.5 around current box
  % keep among neighbors only those within interval
  nbr = nbr(abs(dr)/l < 1.5);
end