close all;
clear all;

% Set up parameters
i = 10;
N = 2^i;
tol = 1e-6;
NG = 10;  % number of Chebyshev pts

k = -N/2:N/2-1;
kk = k(:);

x = (0:N-1)/N;
xx = x(:);

func_name = 'fun0';
switch func_name
    case 'funFT'
        fun = @(x,k)funFT(x,k);
    case 'funIFT'
        fun = @(x,k)funIFT(x,k);
    case 'fun0'
        fun = @(x,k)fun_fio_1D(x,k);
end

f = randn(N,1) + sqrt(-1)*randn(N,1);

tic;
[Factor,Rcomp] = IBF_Cheby(fun,xx,kk,NG,tol);
FactorT = toc;

tic;
ATA = BF_adj_apply(BF_apply(Factor,eye(N)));
ApplyT = toc;

Afun = @(i, j)fun_ATA(ATA, i, j);

[x1,x2] = ndgrid((1:N)/N); 
x = [x1(:) x2(:)]'; 
occ = 64;
rank_or_tol = 1e-8;

F = hifie2(Afun,x,occ,rank_or_tol);

err = snorm(N,@(x)(ATA*x - hifie_mv(F,x)),[],[],1);
err = err/snorm(N,mv,[],[],1);
fprintf('hifie_mv err/time: %10.4e \n',err)






