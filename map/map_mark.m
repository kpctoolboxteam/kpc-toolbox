function MMAP = map_mark(MAP, prob)
% MAP_MARK. Marks arrivals from a MAP according to give probabilities.
% MMAP = MAP_MARK(MAP, PROB) takes a MAP with a single arrival type and
% returns a MMAP with same unmarked inter-arrival process but marked
% arrivals specified by PROB. PROB(k) is the probability that the arrival
% will be marked with type k.

if (sum(prob) < 1 - 1e-6) || (sum(prob) > 1 + 1e-6)
    warning('Input marking probabilities do not sum to 1. Normalizing.');
    prob = prob / sum(prob);
end

MMAP = MAP;
for k=1:length(prob)
    MMAP{2+k} = prob(k) * MAP{2};
end
MMAP = mmap_normalize(MMAP);
end