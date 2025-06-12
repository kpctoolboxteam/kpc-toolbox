function skew = trace_skew(S)
% [skew] = trace_skew(S)
%
% DESCRIPTION
% Compute the bias-corrected skewness of trace S (equivalent to skewness(S,0))

% Remove NaNs
S = S(:);
S = S(~isnan(S));

n = numel(S);
if n < 3
    skew = NaN;
    return;
end

res = S - mean(S, 'omitnan');
s2 = mean(res.^2, 'omitnan');
m3 = mean(res.^3, 'omitnan');

% Bias correction
skew = m3 * (sqrt((n - 1) / n) * n / (n - 2))  / s2^(3/2);
end
