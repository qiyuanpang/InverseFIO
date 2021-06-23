close all;
clear all;

startup;

tol = 1e-12
kernel = "1"
fun_name = "fun" + kernel;
func_name = fun_name;
occ = 32;
rank_or_tol = 1e-3
tol_sol = 1e-8
maxit = 200
repeat_num = 5;
n0 = 8;
tt = 6;
rand_or_cheb = 'cheb';

% dims = [16 25 36]
dims = [16 36 64 100 144 196 256]
cases = length(dims);
apptime = zeros(cases, 1);
bferr = zeros(cases, 1);
condAs = zeros(cases, 1);
condATAs = zeros(cases, 1);

soltime = zeros(cases, 1);
apperr = zeros(cases, 1);
solerr = zeros(cases, 1);

ranks_hqr = zeros(cases, 1);
ranks_hodlr = zeros(cases, 1);
solerrpcg = zeros(cases, 1);
iters = zeros(cases, 1);
for i = 1:cases
    n = dims(i);
    N = n^2;
    NG = 5*floor(log2(N));
    rk = 9*floor(log2(N));
    k = -80/2:80/2-1;
    kk = k(:);
    %rank_or_tol = 8*floor(log2(N))    

    x = (0:N-1)/N;
    xx = x(:);
    

    fprintf('N = %4d \n', N)

    switch fun_name
        case "fun1"
            fun = @(x,k)fun_fio_1D(x,k);
        case "fun2"
            fun = @(x,k)fun_fio2_1D(x,k);
        case "fun3"
            fun = @(x,k)fun_fio3_1D(x,k,0.1);
        case "fun4"
            fun = @(x,k)fun_fio4_1D(N, x, k, 0.05);
        case "fun5"
            fun = @(x,k)fun_fio5_1D(x,k);
        case "fun6"
            fun = @(x,k)fun_fio6_1D(N, x, k, 0.1);
        case "fun7"
            fun = @(x,k)fun_fio6_1D(N, x, k, 0.05);

    end

    A = fun(xx, kk);
    tic;
    for j = 1:repeat_num
        [Y T R] = householderqr(A);
    end
    apptime(i) = toc/repeat_num;
    size(Y)
    size(T)
    size(R)
end

N = dims.^2;
logN = log2(N);
NlogN = logN + log2(logN);
N2logN = logN + 2*log2(logN);

fig = figure(1);
hold on;
h(1) = plot(logN, log2(apptime));
h(2) = plot(logN, logN-logN(1)+log2(apptime(1)));
h(3) = plot(logN, NlogN-NlogN(1)+log2(apptime(1)));
xlabel('Log(N)');
ylabel('Log(Time)/s'); 
title('Forward Time Scaling');
legend(h, 'App time', 'N', 'N log N');
axis square;
saveas(fig, "qr_wybased.png");
hold off;
