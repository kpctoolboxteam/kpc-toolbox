function Q=ctmc_makeinfgen(Q)
A=Q-diag(diag(Q)); Q=A-diag(sum(A,2));
end