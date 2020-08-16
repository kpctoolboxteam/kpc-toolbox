function Prev=dtmc_timereverse(P)
K=length(P);
Prev=P;
pie=dtmc_solve(P);
for i=1:K
    for j=1:K
    Prev(i,j)=P(i,j)*pie(i)/pie(j);
    end
end
Prev=Prev';
end