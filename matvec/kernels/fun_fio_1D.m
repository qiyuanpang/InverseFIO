function res = fun_fio_1D(x,k)

xk = x*k';
sx = (2 + sin(2*pi*x))/8;
tmp = (2*pi)* (xk + sx*abs(k'));
% tmp = (2*pi) * xk;
res = complex(cos(tmp),sin(tmp));

end
