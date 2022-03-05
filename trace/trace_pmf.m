function [pmf,px] = trace_pmf(X)
[pmf] = hist(X, max(X))' ./ numel(X); 
px = unique(X);
end
