close all;
clear all;

startup;

tol = 1e-10
func_name = 'fun1'
occ = 64;
rank_or_tol = 1e-10
repeat_num = 5;

dims = [4:7]
cases = length(dims);
apptime = zeros(cases, 1);
soltime = zeros(cases, 1);
apperr = zeros(cases, 1);
solerr = zeros(cases, 1);
condAs = zeros(cases, 1);
condATAs = zeros(cases, 1);
for i = 1:cases
    ii = dims(i);
    N = 2^(2*ii);
    n = 2^ii;
    NG = 2*ii;
    k = -N/2:N/2-1;
    kk = k(:);

    x = (0:N-1)/N;
    xx = x(:);
    
    fprintf('N = %4d \n', N)

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

    tic;
    for j = 1:repeat_num
        [Factor,Rcomp] = IBF_Cheby(fun,xx,kk,NG,tol);
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

    Afun = @(i, j)fun_ATA(ATA, i, j);

    [x1,x2] = ndgrid((1:n)/n); 
    x = [x1(:) x2(:)]'; 

    F = hifie2my(Afun,x,occ,rank_or_tol);

    err = snorm(N,@(x)(ATA*x - hifie_mv(F,x)),[],[],1);
    err = err/snorm(N,@(x)(ATA*x),[],[],1);
    fprintf('hifie_mv err: %10.4e \n',err)
    apperr(i) = err;

    sol = rand(N,1);
    b = BF_adj_apply(Factor, BF_apply(Factor,sol));
    tic;
    for j = 1:repeat_num
        app = hifie_sv(F, b);
    end
    sol_time = toc/repeat_num;
    err_sol = norm(sol-app)/norm(sol);
    fprintf('hifie_sv err/time: %10.4e/%10.4e (s) \n', err_sol, sol_time)
    solerr(i) = err_sol;
    soltime(i) = sol_time;
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
saveas(fig, "apptime_" + func_name + ".png");

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
saveas(fig, "soltime_" + func_name + ".png");

fig = figure(3);
hold on;
h(1) = plot(logN, log10(apperr));
xlabel('Log(N)');
ylabel('Log10(error)'); 
title('Error of Application');
hold off;
% legend(h, 'App time', 'N log N', 'N log^2 N');
saveas(fig, "apperr_" + func_name + ".png");

fig = figure(4);
hold on;
h(1) = plot(logN, log10(solerr));
xlabel('Log(N)');
ylabel('Log10(error)'); 
title('Error of Solution');
hold off;
% legend(h, 'App time', 'N log N', 'N log^2 N');
saveas(fig, "solerr_"+ func_name + ".png");

fig = figure(5);
h = zeros(1);
h(1) = plot(logN, log10(condAs)); hold on;
xlabel('Log(N)');
ylabel('Log10(cond num)'); 
title('Condition numbers');
hold off;
% legend(h, 'App time', 'N log N', 'N log^2 N');
saveas(fig, "conda_"+ func_name + ".png");

fig = figure(6);
hold on;
h(1) = plot(logN, log10(condATAs));
xlabel('Log(N)');
ylabel('Log10(cond num)'); 
title('Condition numbers');
hold off;
% legend(h, 'App time', 'N log N', 'N log^2 N');
saveas(fig, "condata_" + func_name + ".png");

exit
