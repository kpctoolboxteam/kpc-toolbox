function res=dtmc_isfeasible(P)
sP = sum(P,2);
res = 0;
for tol=1:15
    if min(sP) > 1-10^-tol && max(sP) < 1+10^-tol && min(P(:)) > -10^-tol
        res = tol;
    end
end

end