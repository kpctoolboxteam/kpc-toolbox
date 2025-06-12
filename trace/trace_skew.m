function skew = trace_skewness(S)
% [skew] = trace_skewness(S)
%
% DESCRIPTION
% Compute the bias-corrected skewness of trace S (equivalent to skewness(S,0))

n = numel(S);
m = mean(S);
s = std(S, 1); % population standard deviation

skew = (n / ((n - 1) * (n - 2))) * sum(((S - m) / s).^3);
end
