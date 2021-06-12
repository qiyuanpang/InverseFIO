close all;
clear all;

startup;

tol = 1e-14
kernel = "7"
fun_name = "fun" + kernel;
func_name = fun_name;
occ = 32;
rank_or_tol = 1e-3
tol_sol = 1e-8
maxit = 50
repeat_num = 5;
n0 = 8;
tt = 5;
rand_or_cheb = 'rand';

%dims = [16 25 36]
dims = [16 25 36 49 64 81]
cases = length(dims);
apptime = zeros(cases, 1);
bferr = zeros(cases, 1);
condAs = zeros(cases, 1);
condATAs = zeros(cases, 1);

soltime_hif = zeros(cases, 1);
soltime_hqr = zeros(cases, 1);
apperr_hif = zeros(cases, 1);
apperr_hqr = zeros(cases, 1);
solerr_hif = zeros(cases, 1);
solerr_hqr = zeros(cases, 1);

ranks_hif = zeros(cases, 1);
ranks_hqr = zeros(cases, 1);
solerrpcg_hif = zeros(cases, 1);
solerrpcg_hqr = zeros(cases, 1);
iters_hif = zeros(cases, 1);
iters_hqr = zeros(cases, 1);
for i = 1:cases
    n = dims(i);
    N = n^2;
    NG = 5*floor(log2(N));
    rk = 8*floor(log2(N));
    k = -N/2:N/2-1;
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
            fun = @(x,k)fun_fio3_1D(x,k,1);
        case "fun4"
            fun = @(x,k)fun_fio4_1D(N, x, k, 0.05);
        case "fun5"
            fun = @(x,k)fun_fio5_1D(x,k);
        case "fun6"
            fun = @(x,k)fun_fio6_1D(N, x, k, 0.1);
        case "fun7"
            fun = @(x,k)fun_fio6_1D(N, x, k, 0.05);

    end

    tic;
    for j = 1:repeat_num
        %[Factor,Rcomp] = IBF_Cheby(fun,xx,kk,NG,tol);
        Factor = BF_IDBF(fun, xx, kk, n0, rk, tol, rand_or_cheb, tt);
    end
    FactorT = toc/repeat_num;

    apptime(i) = FactorT;

    tic;
    A = BF_apply(Factor,eye(N));
    ATA = BF_adj_apply(Factor, A);
    ApplyT = toc;

    condA = cond(A);
    condATA = cond(ATA);
    fprintf('condition number estimation of A  : %10.4e \n', condA)
    fprintf('condition number estimation of ATA: %10.4e \n', condATA)

    condAs(i) = condA;
    condATAs(i) = condATA;
    
    A1 = fun(xx, kk);
    ATA1 = A1'*A1;

    condA1 = cond(A1);
    condATA1 = cond(ATA1);
    fprintf('condition number estimation of A1  : %10.4e \n', condA1)
    fprintf('condition number estimation of ATA1: %10.4e \n', condATA1)
    
    randx = rand(N, 1);
    bfacc = norm(A*randx-A1*randx)/norm(A1*randx);
    fprintf('err between matrices: %10.4e \n', bfacc)
    %fprintf('err between ATAs    : %10.4e \n', norm(ATA*randx-ATA1*randx)/norm(ATA1*randx))
    bferr(i) = bfacc;    

    Afun = @(i, j)fun_ATA(ATA, i, j);

    [x1,x2] = ndgrid((1:n)/n); 
    x = [x1(:) x2(:)]'; 

    [F, rk] = hifie2my(Afun,x,occ,rank_or_tol);

    fprintf('rank of hif: %6d \n', rk)
    ranks_hif(i) = rk;

    [Y, T, R, rk] = hqr(ATA, [], [], [], floor(log2(N)/2)+1, rank_or_tol);

    fprintf('rank of hqr: %6d \n', rk)
    ranks_hqr(i) = rk;

    err = snorm(N,@(x)(ATA*x - hifie_mv(F,x)),[],[],1,128);
    err = err/snorm(N,@(x)(ATA*x),[],[],1, 128);
    fprintf('hifie_mv err: %10.4e \n',err)
    apperr_hif(i) = err;
    
    Q = eye(N)-Y*T*Y';
    R = [R;zeros(N-size(R,1), size(R,2))];
    x = rand(N,1);
    err_app = norm(Q*R*x-ATA*x)/norm(ATA*x);
    fprintf('hqr_mv   err: %10.4e \n', err_app)
    apperr_hqr(i) = err_app;
    
    sol = 2*rand(N,1);
    b = BF_adj_apply(Factor, A1*sol);
    %b = ATA1*sol;
    tic;
    for j = 1:repeat_num
        app = hifie_sv(F, b);
    end
    sol_time = toc/repeat_num;
    err_sol = norm(sol-app)/norm(sol);
    fprintf('hifie_sv err/time: %10.4e/%10.4e (s) \n', err_sol, sol_time)
    solerr_hif(i) = err_sol;
    soltime_hif(i) = sol_time;
    
    tic;
    for j = 1:repeat_num
        app = R\(b - Y*T'*Y'*b);
    end
    sol_time = toc/repeat_num;
    err_sol = norm(sol-app)/norm(sol);
    fprintf('hqr_sv   err/time: %10.4e/%10.4e (s) \n', err_sol, sol_time)
    solerr_hqr(i) = err_sol;
    soltime_hqr(i) = sol_time;
    
    b = rand(N, 1);
    pcd = @(x)hifie_sv(F, x);
    [x, flag, relres, iter] = pcg(ATA, b, tol_sol, maxit, pcd);
    solerrpcg_hif(i) = relres;
    iters_hif(i) = iter;
    fprintf('Solve the equation by PCG with HIF as a preconditioner in %4d iterations, rel error: %10.4e \n', iter, relres)
    
    pcd = @(b)R\(b - Y*T'*Y'*b);
    [x, flag, relres, iter] = pcg(ATA, b, tol_sol, maxit, pcd);
    solerrpcg_hqr(i) = relres;
    iters_hqr(i) = iter;
    fprintf('Solve the equation by PCG with HQR as a preconditioner in %4d iterations, rel error: %10.4e \n', iter, relres)
    
end

N = dims.^2;
logN = log2(N);
NlogN = logN + log2(logN);
N2logN = logN + 2*log2(logN);


fig = figure(1);
hold on;
h(1) = plot(logN, log2(apptime));
h(2) = plot(logN, NlogN-NlogN(1)+log2(apptime(1)));
h(3) = plot(logN, N2logN-N2logN(1)+log2(apptime(1)));
xlabel('Log(N)');
ylabel('Log(Time)/s'); 
title('BF time scaling');
legend(h, 'App time', 'N log N', 'N log^2 N');
axis square;
saveas(fig, "./comp/1D/kernel" + kernel + "/apptime_" + func_name + ".png");
hold off;

fig = figure(2);
hold on;
h(1) = plot(logN, log2(soltime_hif));
h(2) = plot(logN, log2(soltime_hqr));
h(3) = plot(logN, NlogN-NlogN(1)+(log2(soltime_hif(1))+log2(soltime_hqr(1)))/2);
h(4) = plot(logN, N2logN-N2logN(1)+(log2(soltime_hif(1))+log2(soltime_hqr(1)))/2);
xlabel('Log(N)');
ylabel('Log(Time)/s'); 
title('Equation Solving time scaling');
legend(h, 'HIF', 'HQR', 'N log N', 'N log^2 N');
saveas(fig, "./comp/1D/kernel" + kernel + "/soltime_" + func_name + ".png");
hold off;

fig = figure(3);
hold on;
h(1) = plot(logN, log10(apperr_hif));
h(2) = plot(logN, log10(apperr_hqr));
xlabel('Log(N)')
ylabel('Log10(error)'); 
title('Forward error');
legend(h, 'HIF', 'HQR');
saveas(fig, "./comp/1D/kernel" + kernel + "/apperr_" + func_name + ".png");
hold off;

fig = figure(4);
hold on;
h(1) = plot(logN, log10(solerr_hif));
h(2) = plot(logN, log10(solerr_hqr));
xlabel('Log(N)');
ylabel('Log10(error)'); 
title('Backward error');
legend(h, 'HIF', 'HQR');
saveas(fig, "./comp/1D/kernel" + kernel + "/solerr_"+ func_name + ".png");
hold off;

fig = figure(5);
hold on;
h(1) = plot(logN, log10(condAs));
h(2) = plot(logN, log10(condATAs));
xlabel('Log(N)');
ylabel('Log10(cond num)'); 
title('Condition numbers');
legend(h, 'A', 'ATA');
saveas(fig, "./comp/1D/kernel" + kernel + "/conda_"+ func_name + ".png");
hold off;

fig = figure(6);
hold on;
h(1) = plot(logN, log10(bferr));
xlabel('Log(N)');
ylabel('Log10(error)'); 
title('Error of BF');
% legend(h, 'App time', 'N log N', 'N log^2 N');
saveas(fig, "./comp/1D/kernel" + kernel + "/bferr_"+ func_name + ".png");
hold off;

fig = figure(7);
hold on;
h(1) = plot(logN, log2(ranks_hif));
h(2) = plot(logN, log2(ranks_hqr));
h(3) = plot(logN, 0.5*logN-0.5*logN(1)+(log2(ranks_hif(1))+log2(ranks_hqr(1)))/2);
h(4) = plot(logN, log2(logN)-log2(logN(1))+(log2(ranks_hif(1))+log2(ranks_hqr(1)))/2);
xlabel('Log(N)');
ylabel('Log2(rank)'); 
title('Rank');
legend(h, 'HIF', 'HQR', 'sqrt(N)', 'log N');
saveas(fig, "./comp/1D/kernel" + kernel + "/rank_"+ func_name + ".png");
hold off;

fig = figure(8);
hold on;
h(1) = plot(logN, log10(solerrpcg_hif));
h(2) = plot(logN, log10(solerrpcg_hqr));
xlabel('Log(N)');
ylabel('Log10(error)'); 
title('Error of PCG Solution');
legend(h, 'HIF', 'HQR');
saveas(fig, "./comp/1D/kernel" + kernel + "/pcgerr_"+ func_name + ".png");
hold off;

fig = figure(9);
hold on;
h(1) = plot(logN, iters_hif);
h(2) = plot(logN, iters_hqr);
xlabel('Log(N)');
ylabel('Iteration'); 
title('Iterations in PCG');
legend(h, 'HIF', 'HQR');
saveas(fig, "./comp/1D/kernel" + kernel + "/iters_" + func_name + ".png");
hold off;

exit;
