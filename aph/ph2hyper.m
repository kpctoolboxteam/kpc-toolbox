function [lambda,prob]=ph2hyper(PH)
if (norm(PH{1}-diag(diag(PH{1})))>1e-10)
    error('The PH distribution is not hyper-exponential')
end
lambda=diag(-PH{1})';
prob=dtmc_solve(inv(-PH{1})*PH{2});
end

