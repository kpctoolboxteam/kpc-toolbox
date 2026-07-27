function [sts]=dtmc_simulate(P, pi0, n)
% [sts]=dtmc_simulate(P, pi0, n)

rnd  = rand;
cpi0  = cumsum(pi0);
if isempty(intersect(find( rnd - cpi0 > 0), find(pi0>0)))
    st = min(find(pi0));
else
    st = min(intersect(1 + find( rnd - cpi0 > 0),find(pi0>0)));
end

F = cumsum(P,2);
for i=1:n
    sts(i) = st; soujt(i)=1;
    if F(st,end)==0 || P(st,st)==1
        return;
    end
    rnd = rand;
    if isempty(intersect(find( rnd - F(st,:) > 0),find(P(st,:)>0)))
        st = min(find(P(st,:)>0));
    else
        st = min(intersect(1 + find( rnd - F(st,:) > 0),find(P(st,:)>0)));
    end
end
end