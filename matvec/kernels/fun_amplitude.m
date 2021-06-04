function ans = fun_amplitude(x, k, xk, kk, sigma2)
nk = length(xk);
n = length(x);
m = length(k);
matx = zeros(n,nk);
matk = zeros(m,nk);
%fprintf('size %d %d \n', length(xk), length(kk))
for i  = 1:nk
   %fprintf('%d %10.4e %10.4e \n',i, xk(i), kk(i))
   matx(:, i) = (x - xk(i)).^2;
   matk(:, i) = (k - kk(i)).^2;
end
ans = zeros(n,m);
for i = 1:n
   for j = 1:m
       ans(i,j) = sum(exp(-(matx(i,:)+matk(j,:))./sigma2));
   end
end
end
