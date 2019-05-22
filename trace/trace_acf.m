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

if exist('xcorr')
    rho=xcorr(S-mean(S),max(lags),'coeff');
    rho=rho((max(lags)+2):end);
    rho=rho(lags);
else
    %warning('trace_acf: signal processing toolbox not found, autocorrelation estimation may be slower than usual.');
    S=S(:);
    acv = autocov(S-mean(S));
    rho=acv(1+lags)'/acv(1);
%     E1=mean(S);
%     E2=mean(S.^2);
%     rho = zeros(1,length(lags));
%     for ki=1:length(lags)
%         rho(ki)=(trace_joint(S,[0,lags(ki)],[1,1])-E1^2)/(E2-E1^2);
%     end
end

rho=rho(:);

end