function rho=trace_acf(S,lags)
% [rho]=trace_acf(S,L)
%
% DESCRIPTION
% Compute the autocorrelation function for trace S at the lags specified in vector L

if nargin==1
    lags=1;
end

rho=[];

if max(lags)>length(S)-2
    warning('Lags are higher than the trace length, truncating the lags vector.');
    lags = min(lags,length(S)-2);
    lags(lags == length(S)-2) = [];
end

% Always estimate through autocov(), the unbiased 1/(N-k) estimator. The
% Signal Processing Toolbox xcorr(...,'coeff') path that used to be taken
% whenever xcorr was on the path is the biased 1/N estimator, so the two
% disagreed by O(k/N) and the answer depended on which toolboxes happened to
% be installed.
S=S(:);
acv = autocov(S-mean(S));
rho=acv(1+lags)'/acv(1);

rho=rho(:);

end