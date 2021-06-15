function res = fun0var(N,x,k,s)

xk = x*k';
sx = (2 + sin(2*pi*x))/7;
tmp = (2*pi)* (xk + sx*abs(k'));
[K,X] = meshgrid(k,x);

a = exp( -((X-0.5).^2+(K/N).^2)/s)+exp( -((X-0.6).^2+(K/N-0.4).^2)/s)+exp( -((X-0.9).^2+(K/N-0.3).^2)/s)+exp( -((X+0.2).^2+(K/N+0.3).^2)/s)...
    +exp( -((X-0.5-1).^2+(K/N).^2)/s)+exp( -((X-0.6-1).^2+(K/N-0.4).^2)/s)+exp( -((X-0.9-1).^2+(K/N-0.3).^2)/s)+exp( -((X+0.2-1).^2+(K/N+0.3).^2)/s)...
    +exp( -((X-0.5+1).^2+(K/N).^2)/s)+exp( -((X-0.6+1).^2+(K/N-0.4).^2)/s)+exp( -((X-0.9+1).^2+(K/N-0.3).^2)/s)+exp( -((X+0.2+1).^2+(K/N+0.3).^2)/s);

res = a.*complex(cos(tmp),sin(tmp));

end
