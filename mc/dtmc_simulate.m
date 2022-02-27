function [sts]=dtmc_simulate(P, pi0, n)
% [sts]=dtmc_simulate(P, pi0, n)

[~,st] =  min(abs(rand-cumsum(pi0)));
F = cumsum(P,2);
for i=1:n
    sts(i) = st; soujt(i) =1;
    st =  1 + max([0,find( rand - F(st,:) > 0)]);
end

end