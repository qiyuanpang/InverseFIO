function res = fun_fio5_1D(x,k)

xk = x*k';
sx = (2 + sin(2*pi*x))/7;
tmp = (2*pi)* (xk + sx*abs(k)');
res = complex(cos(tmp),sin(tmp));

end
