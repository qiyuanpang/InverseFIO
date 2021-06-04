function res = fun_fio6_1D(x,k,ax,ak,sigma2)

xk = x*k';
sx = (2 + sin(2*pi*x))/8;
tmp = (2*pi)* (xk + sx*abs(k'));
res = complex(cos(tmp),sin(tmp));

end
