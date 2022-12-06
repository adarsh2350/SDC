function f = fun(x,lb,ub,ht)
f = -0.001*sum(ht.*(x - lb)./((1+ht)./(ub - x)));
