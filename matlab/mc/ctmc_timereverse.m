function Qrev=ctmc_timereverse(Q)
K=length(Q);
Qrev=Q;
pie=ctmc_solve(Q);
for i=1:K
    for j=1:K
    Qrev(i,j)=Q(i,j)*pie(i)/pie(j);
    end
end
Qrev=Qrev';
end