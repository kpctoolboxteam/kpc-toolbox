function APH=aph_rand(K)
e = ones(K,1);
P=e*dtmc_solve(dtmc_rand(K));
mu=rand(1,K);
D0=-diag(mu)+diag(mu(1:K-1),1);
D1=-D0*P;
APH={D0,D1};
end

