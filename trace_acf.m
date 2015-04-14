function rho=trace_acf(S,lags)
% [rho]=trace_acf(S,L)
%
% DESCRIPTION
% Compute the autocorrelation function for trace S at the lags specified in vector L 

if nargin==1
    lags=1;
end


rho=[];

if exist('xcorr')
    rho=xcorr(S-mean(S),max(lags),'coeff');
    rho=rho((max(lags)+2):end);
    rho=rho(lags);
else
    warning('trace_acf: signal processing toolbox not found, autocorrelation estimation may be slower than usual.');
    E1=mean(S);
    E2=mean(S.^2);
    for k=lags
        rho(end+1)=(trace_joint(S,[0,k],[1,1])-E1^2)/(E2-E1^2);
    end
end

rho=rho(:);

end