function P=ctmc_randomization(Q,q)
if nargin==1
    q=(max(max(abs(Q))))+rand;
end
P=Q/q + eye(size(Q));
end