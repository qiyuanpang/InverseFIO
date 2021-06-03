close all;
clear all;

startup;

tol = 1e-15
func_name = 'fun1'
occ = 32;
%rank_or_tol = 1e-6
tol_sol = 1e-8
maxit = 50
repeat_num = 5;
n0 = 8;
tt = 5;
rand_or_cheb = 'rand';

dims = [4:7]
cases = length(dims);
apptime = zeros(cases, 1);
soltime = zeros(cases, 1);
apperr = zeros(cases, 1);
solerr = zeros(cases, 1);
condAs = zeros(cases, 1);
condATAs = zeros(cases, 1);
bferr = zeros(cases, 1);
ranks = zeros(cases, 1);
solerrpcg = zeros(cases, 1);
iters = zeros(cases, 1);
for i = 1:cases
    ii = dims(i);
    N = 2^(2*ii);
    n = 2^ii;
    NG = 5*ii;
    rk = 15*ii;
    rank_or_tol = 15*ii
    k = -n/2:n/2-1;
    [k1,k2] = ndgrid(k);
    kk = [k1(:) k2(:)];

    x = (0:n-1)/n;
    [x1,x2] = ndgrid(x);
    xx = [x1(:) x2(:)];
        
    fprintf('N = %4d \n', N)

    switch func_name
        case 'fun1'
            fun = @(x,k)fun_fio_2D(x,k);
        case 'fun2'
            fun = @(x,k)fun_fio2_2D(x,k);
        case 'fun3'
            fun = @(x,k)fun_fio3_2D(x,k);
        case 'fun4'
            fun = @(x,k)fun_fio4_2D(x,k);

    end

    tic;
    for j = 1:repeat_num
        [Factor,Rcomp] = IBF_Cheby(fun,xx,kk,NG,tol);
        % Factor = BF_IDBF(fun, xx, kk, n0, rk, tol, rand_or_cheb, tt);
    end
    FactorT = toc/repeat_num;

    apptime(i) = FactorT;

    tic;
    A = BF_apply(Factor,eye(N));
    ATA = BF_adj_apply(Factor, A);
    ApplyT = toc;

    condA = condest(A);
    condATA = condest(ATA);
    fprintf('condition number estimation of A  : %10.4e \n', condA)
    fprintf('condition number estimation of ATA: %10.4e \n', condATA)

    condAs(i) = condA;
    condATAs(i) = condATA;
    
    A1 = fun(xx, kk);
    ATA1 = A1'*A1;

    condA1 = condest(A1);
    condATA1 = condest(ATA1);
    fprintf('condition number estimation of A1  : %10.4e \n', condA1)
    fprintf('condition number estimation of ATA1: %10.4e \n', condATA1)
    
    randx = rand(N, 1);
    bfacc = norm(A*randx-A1*randx)/norm(A1*randx);
    fprintf('err between matrices: %10.4e \n', bfacc)
    %fprintf('err between ATAs    : %10.4e \n', norm(ATA*randx-ATA1*randx)/norm(ATA1*randx))
    bferr(i) = bfacc;    

    Afun = @(i, j)fun_ATA(ATA, i, j);

    % [x1,x2] = ndgrid((1:n)/n); 
    % x = [x1(:) x2(:)]'; 

    [F, rank] = hifie2my(Afun,xx',occ,rank_or_tol);

    fprintf('rank of hif: %6d \n', rank)
    ranks(i) = rank;

    err = snorm(N,@(x)(ATA*x - hifie_mv(F,x)),[],[],1,128);
    err = err/snorm(N,@(x)(ATA*x),[],[],1, 128);
    fprintf('hifie_mv err: %10.4e \n',err)
    apperr(i) = err;

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
    solerr(i) = err_sol;
    soltime(i) = sol_time;

    %function y = pcd(x)
    %   y = hifie_sv(F, BF_adj_apply(Factor, x));
    %end
    b = rand(N, 1);
    pcd = @(x)hifie_sv(F, x);
    [x, flag, relres, iter] = pcg(ATA, b, tol_sol, maxit, pcd);
    solerrpcg(i) = relres;
    iters(i) = iter;
    fprintf('Solve the equation by PCG with F as a preconditioner in %4d iterations, rel error: %10.4e \n', iter, relres)
    
    %pcd = @(x)hifie_sv(F, x);
    [x, flag, relres, iter] = pcg(ATA, b, tol_sol, maxit);
    fprintf('Solve the equation by PCG without preconditioners    in %4d iterations, rel error: %10.4e \n', iter, relres)

    %[x, flag, relres, iter] = pcg(A, b, tol_sol, 100);
    %fprintf('Solve the original equation by PCG                   in %4d iterations, rel error: %10.4e \n', iter, relres)
end

N = 2.^(2*dims);
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
title('Time scaling of Application');
legend(h, 'App time', 'N log N', 'N log^2 N');
axis square;
saveas(fig, "apptime_" + func_name + "_2D.png");

fig = figure(2);
hold on;
h(1) = plot(logN, log2(soltime));
h(2) = plot(logN, NlogN-NlogN(1)+log2(soltime(1)));
h(3) = plot(logN, N2logN-N2logN(1)+log2(soltime(1)));
xlabel('Log(N)');
ylabel('Log(Time)/s'); 
title('Time scaling of Solving equations');
hold off;
legend(h, 'Sol time', 'N log N', 'N log^2 N');
saveas(fig, "soltime_" + func_name + "_2D.png");

fig = figure(3);
hold on;
h(1) = plot(logN, log10(apperr));
xlabel('Log(N)');
ylabel('Log10(error)'); 
title('Error of Application');
hold off;
% legend(h, 'App time', 'N log N', 'N log^2 N');
saveas(fig, "apperr_" + func_name + "_2D.png");

fig = figure(4);
hold on;
h(1) = plot(logN, log10(solerr));
xlabel('Log(N)');
ylabel('Log10(error)'); 
title('Error of Solution');
hold off;
% legend(h, 'App time', 'N log N', 'N log^2 N');
saveas(fig, "solerr_"+ func_name + "_2D.png");

fig = figure(5);
hold on;
h(1) = plot(logN, log10(condAs));
xlabel('Log(N)');
ylabel('Log10(cond num)'); 
title('Condition numbers');
hold off;
% legend(h, 'App time', 'N log N', 'N log^2 N');
saveas(fig, "conda_"+ func_name + "_2D.png");

fig = figure(6);
hold on;
h(1) = plot(logN, log10(condATAs));
xlabel('Log(N)');
ylabel('Log10(cond num)'); 
title('Condition numbers');
hold off;
% legend(h, 'App time', 'N log N', 'N log^2 N');
saveas(fig, "condata_" + func_name + "_2D.png");

fig = figure(7);
hold on;
h(1) = plot(logN, log10(bferr));
xlabel('Log(N)');
ylabel('Log10(error)'); 
title('Error of BF');
hold off;
% legend(h, 'App time', 'N log N', 'N log^2 N');
saveas(fig, "bferr_"+ func_name + "_2D.png");


fig = figure(8);
hold on;
h(1) = plot(logN, log2(ranks));
h(2) = plot(logN, 0.5*logN-0.5*logN(1)+log2(ranks(1)));
h(3) = plot(logN, log2(logN)-log2(logN(1))+log2(ranks(1)));
xlabel('Log(N)');
ylabel('Log2(rank)'); 
title('Rank');
hold off;
legend(h, 'rank', 'sqrt(N)', 'log N');
saveas(fig, "rank_"+ func_name + "_2D.png");

fig = figure(9);
hold on;
h(1) = plot(logN, log10(solerrpcg));
xlabel('Log(N)');
ylabel('Log10(error)'); 
title('Error of Solution via PCG');
hold off;
% legend(h, 'App time', 'N log N', 'N log^2 N');
saveas(fig, "pcgerr_"+ func_name + "_2D.png");

fig = figure(10);
hold on;
h(1) = plot(logN, iters);
xlabel('Log(N)');
ylabel('Iteration'); 
title('Iterations in PCG');
hold off;
% legend(h, 'App time', 'N log N', 'N log^2 N');
saveas(fig, "iters_" + func_name + "_2D.png");


exit
