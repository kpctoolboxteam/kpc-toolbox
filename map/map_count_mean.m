function [m] = map_count_mean(map,t)
% Computes the mean of the counting process, at resolution t, for the
% given MAP.
% Input:
% - mmap: Markovian Arrival Process
% - t: the period considered for each sample of the counting process
% Output:
% - m: mean arrivals in (0,t]

n = size(map{1},1);

e = ones(n,1);
theta = map_prob(map);

% arrival rate
l = theta * map{2} * e;

% mean
if (size(t,1) > 1)
    m = l * t;
else
    m = l * t';
end

end