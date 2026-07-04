clear;

H{1} = aph_from2moments(1,1.3);
K=7;

pi0 = map_pie(H{1});

for j=2:floor(K/2)
    H{j} = map_hyperexp(2^(j-1),1.3^(j-1),0.99);
    pi0 = kron(pi0,map_pie(H{j}));
end

for j=floor(K/2)+1:K
    H{j} = map_hyperexp(2^(j-1),1.3^(j-1),0.99);
    pi0 = kron(pi0,map_pie(H{j}));
end

S = map_kpc(H);
map_isfeasible(S)
target_scv=map_scv(S)

kpcsample = [];
while isempty(kpcsample) || norm(var(kpcsample)/mean(kpcsample).^2-target_scv,1)>0.10*target_scv
    kpcsample = [map_sample(S,1e6,pi0,2300)];
    empirical_scv = var(kpcsample)/mean(kpcsample).^2
end
empirical_scv = var(kpcsample)/mean(kpcsample).^2

kpcsample=kpcsample/mean(kpcsample);
