function f = fun(x,lb,ub,ht)
f = sum(((ub - x).*ht) + ((ht/75).*(x - lb)));
