function Sred = map_shrink(S,newsize)

pi = ctmc_solve(S{1} + S{2});
alphav = map_pie(S);

[~,indL] = maxk(pi,newsize);
Aind = sort(indL);

alphacompl = alphav(Aind)./sum(alphav(Aind));
Sred{1} = S{1}(Aind,Aind);
Sred{2} = -sum(Sred{1},2)*alphacompl;

end