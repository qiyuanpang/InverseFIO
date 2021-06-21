startup;

n = 25

n0 = 8;
tt = 6;
rand_or_cheb = 'cheb';
N = n^2;
NG = 5*floor(log2(N));
rk = 9*floor(log2(N));
k = -N/2:N/2-1;
kk = k(:);
tol = 1e-8;
tol_or_rank = 1E-6; 

x = (0:N-1)/N;
xx = x(:);
fun = @(x,k)fun_fio_1D(x, k);

Factor = BF_IDBF(fun, xx, kk, n0, rk, tol, rand_or_cheb, tt);

A = BF_apply(Factor,eye(N));
ATA = BF_adj_apply(Factor, A);


lvls = 4;
F = hodlr(ATA, lvls, tol);
norm2 = norm(ATA);
[YA, YB, YC, T, R] = hodlrqr(F, [], [], [], lvls, norm2, tol_or_rank);
Y = YA;
v = rand(N, 1);
y = hodlr_apply(R, v);
y = y - hodlr_apply(Y, hodlr_apply(T, hodlr_adj_apply(Y, y)));
norm(y - ATA*v)/norm(ATA*v)
