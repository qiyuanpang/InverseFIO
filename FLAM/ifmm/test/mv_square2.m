% Source-target evaluations on the unit square, Helmholtz kernel.
%
% This is basically the same as MV_SQUARE1 but using the Helmholtz kernel. The
% associated matrix is rectangular and complex.
%
% Inputs (defaults are used if not provided or set empty):
%
%   - M: number of row points (default: M = 16384)
%   - N: number of column points (default: N = 8192)
%   - K: wavenumber (default: K = 2*PI*8)
%   - OCC: tree occupancy parameter (default: OCC = 128)
%   - P: number of proxy points (default: P = 64)
%   - RANK_OR_TOL: local precision parameter (default: RANK_OR_TOL = 1e-9)
%   - TMAX: ID interpolation matrix entry bound (default: TMAX = 2)
%   - NEAR: near-field compression parameter (default: NEAR = 0)
%   - STORE: storage parameter (default: STORE = 'N')

function mv_square2(M,N,k,occ,p,rank_or_tol,Tmax,near,store)

  % set default parameters
  if nargin < 1 || isempty(M), M = 16384; end
  if nargin < 2 || isempty(N), N =  8192; end
  if nargin < 3 || isempty(k), k = 2*pi*8; end
  if nargin < 4 || isempty(occ), occ = 128; end
  if nargin < 5 || isempty(p), p = 64; end
  if nargin < 6 || isempty(rank_or_tol), rank_or_tol = 1e-9; end
  if nargin < 7 || isempty(Tmax), Tmax = 2; end
  if nargin < 8 || isempty(near), near = 0; end
  if nargin < 9 || isempty(store), store = 'n'; end

  % initialize
  rx = rand(2,M);                                              % row points
  cx = rand(2,N);                                              % col points
  theta = (1:p)*2*pi/p; proxy = 1.5*[cos(theta); sin(theta)];  % proxy points
  % reference proxy points are for unit box [-1, 1]^2

  % compress matrix
  Afun = @(i,j)Afun_(i,j,rx,cx,k);
  pxyfun = @(rc,rx,cx,slf,nbr,l,ctr)pxyfun_(rc,rx,cx,slf,nbr,l,ctr,proxy,k);
  opts = struct('Tmax',Tmax,'near',near,'store',store,'verb',1);
  tic; F = ifmm(Afun,rx,cx,occ,rank_or_tol,pxyfun,opts); t = toc;
  w = whos('F'); mem = w.bytes/1e6;
  fprintf('ifmm time/mem: %10.4e (s) / %6.2f (MB)\n',t,mem)

  % test matrix apply accuracy
  X = rand(N,1); X = X/norm(X);
  tic; ifmm_mv(F,X,Afun,'n'); t = toc;  % for timing
  X = rand(N,16); X = X/norm(X);  % test against 16 vectors for robustness
  r = randperm(M); r = r(1:min(M,128));  % check up to 128 rows in result
  Y = ifmm_mv(F,X,Afun,'n');
  Z = Afun(r,1:N)*X;
  err = norm(Z - Y(r,:))/norm(Z);
  fprintf('ifmm_mv:\n')
  fprintf('  multiply err/time: %10.4e / %10.4e (s)\n',err,t)

  % test matrix adjoint apply accuracy
  X = rand(M,1); X = X/norm(X);
  tic; ifmm_mv(F,X,Afun,'c'); t = toc;  % for timing
  X = rand(M,16); X = X/norm(X);  % test against 16 vectors for robustness
  r = randperm(N); r = r(1:min(N,128));  % check up to 128 rows in result
  Y = ifmm_mv(F,X,Afun,'c');
  Z = Afun(1:M,r)'*X;
  err = norm(Z - Y(r,:))/norm(Z);
  fprintf('  adjoint multiply err/time: %10.4e / %10.4e (s)\n',err,t)
end

% kernel function
function K = Kfun(x,y,k)
  dx = x(1,:)' - y(1,:);
  dy = x(2,:)' - y(2,:);
  K = 0.25i*besselh(0,1,k*sqrt(dx.^2 + dy.^2));
end

% matrix entries
function A = Afun_(i,j,rx,cx,k)
  A = Kfun(rx(:,i),cx(:,j),k);
end

% proxy function
function [Kpxy,nbr] = pxyfun_(rc,rx,cx,slf,nbr,l,ctr,proxy,k)
  pxy = proxy.*l + ctr;  % scale and translate reference points
  if rc == 'r'
    Kpxy = Kfun(rx(:,slf),pxy,k);
    dr = cx(:,nbr) - ctr;
  else
    Kpxy = Kfun(pxy,cx(:,slf),k);
    dr = rx(:,nbr) - ctr;
  end
  % proxy points form ellipse of scaled "radius" 1.5 around current box
  % keep among neighbors only those within ellipse
  nbr = nbr(sum((dr./l).^2) < 1.5^2);
end