function skew = trace_skewness(S)
% [skew] = trace_skewness(S)
%
% DESCRIPTION
% Compute the skewness of trace S

skew = mean((S - mean(S)).^3) / std(S)^3;
end