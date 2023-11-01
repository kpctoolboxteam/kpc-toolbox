function H=hyper_rand(k)
v=rand(k,1);
alpha = rand(1,k); alpha = alpha/sum(alpha);
H0=-diag(v);
H1=diag(v)*ones(k,1)*alpha;
H={H0,H1};
end