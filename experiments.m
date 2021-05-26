
close all;
clear all;

startup;

% Set up parameters
i = 10;
N = 2^i;
n = 2^(i/2);
tol = 1e-6;
NG = 10;  % number of Chebyshev pts

k = -N/2:N/2-1;
kk = k(:);

x = (0:N-1)/N;
xx = x(:);

func_name = 'fun1';
switch func_name
    case 'fun1'
        fun = @(x,k)fun_fio_1D(x,k);
    case 'fun2'
        fun = @(x,k)fun_fio2_1D(x,k);
    case 'fun3'
        fun = @(x,k)fun_fio3_1D(x,k,1);
    case 'fun4'
        fun = @(x,k)fun_fio4_1D(x,k,1);
    case 'fun5'
        fun = @(x,k)fun_fio5_1D(x,k);
end

f = randn(N,1) + sqrt(-1)*randn(N,1);

tic;
[Factor,Rcomp] = IBF_Cheby(fun,xx,kk,NG,tol);
FactorT = toc;

size(BF_apply(Factor, eye(N)));

tic;
A = BF_apply(Factor,eye(N));
ATA = BF_adj_apply(Factor, A);
ApplyT = toc;

fprintf('condition number of A  : %10.4e \n', cond(A))
fprintf('condition number of ATA: %10.4e \n', cond(ATA))

Afun = @(i, j)fun_ATA(ATA, i, j);

[x1,x2] = ndgrid((1:n)/n); 
x = [x1(:) x2(:)]'; 
occ = 64;
rank_or_tol = 1e-8;

F = hifie2my(Afun,x,occ,rank_or_tol);

err = snorm(N,@(x)(ATA*x - hifie_mv(F,x)),[],[],1);
err = err/snorm(N,@(x)(ATA*x),[],[],1);
fprintf('hifie_mv err: %10.4e \n',err)

sol = rand(N,1);
b = BF_adj_apply(Factor, BF_apply(Factor,sol));
tic;
app = hifie_sv(F, b);
sol_time = toc;
err_sol = norm(sol-app)/norm(sol);
fprintf('hifie_sv err/time: %10.4e/%10.4e (s) \n', err_sol, sol_time)






