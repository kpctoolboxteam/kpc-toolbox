function sample = kpcqbd_sample(H, nsamples)
K = length(H);
sample = map_sample(H{1},nsamples);
for j=2:K
    sample = sample .* det_sample(H{j},nsamples);
end
end