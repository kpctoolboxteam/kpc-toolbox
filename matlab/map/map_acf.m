function ACFCOEFFS=map_acf(MAP,LAGS)
% ACFCOEFFS=map_acf(MAP,LAGS) - Compute autocorrelation coefficients of
% interarrival times
%
%  Input:
%  MAP: a MAP in the form of {D0,D1}
%  LAGS: a set of positive integers specificying the lags of the
%  autocorrelation coefficients (DEFAULT: 1)
%
%  Output:
%  ACFCOEFFS: autocorrelation coefficients returned in the same order of
%  the LAGS vector
%
%  Examples:
%  - map_acf(MAP) returns the lag-1 autocorrelation coefficient
%  - map_acf(MAP,1) returns the lag-1 autocorrelation coefficient
%  - map_acf(MAP,1:10) returns the first ten autocorrelation coefficients
%  - map_acf(MAP,logspace(0,4,5)) returns five logarithmically spaced
%  autocorrelation coefficients in [1e0,1e4]
%


if nargin==1
    LAGS=1;
end
P=map_embedded(MAP);
x=map_lambda(MAP)*map_prob(MAP);
if issym(MAP{1}(1,1))    
    if ~isdeployed
        ACFCOEFFS=sym([]);
        y=inv(-MAP{1})*sym(ones(size(MAP{1},1),1));
    end
else
    ACFCOEFFS=([]);
    y=inv(-MAP{1})*ones(size(MAP{1},1),1);
end
for lag=LAGS
    ACFCOEFFS(end+1)=x*(P^lag)*y;
end
ACFCOEFFS=(ACFCOEFFS-1)./map_scv(MAP);
end