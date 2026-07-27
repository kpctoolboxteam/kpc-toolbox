function [GAMMA, RHO0, RESIDUALS] = trace_gamma(T, limit)
% Estimates the auto-correlation decay rate of a trace.
% Input:
% - T: the trace
% - limit: maximum lag considered (optional, default = 1000)
% Output:
% - GAMMA: autocorrelation decay rate

if nargin < 2
    limit = 1000;
end

M1 = mean(T);
M2 = mean(T.^2);

lag = 1:min(limit,length(T));
rho = trace_acf(T, lag);
VAR = M2-M1^2;
SCV = VAR/M1^2;
RHO0 = 1/2 * (1 - 1/SCV);

opt = statset('nlinfit');
opt.MaxIter = 1e5;
opt.Display = 'off';
opt.RobustWgtFun = 'fair';
try
    [GAMMA,RESIDUALS] = nlinfit(lag, rho, @geometric, 0.99, opt);
catch err
    warning('Non linear regression for ACF decay rate failed, trying lsqcurvefit'); 
    GAMMA = lsqcurvefit(@geometric, 0.99, lag, rho, -1, 1);
end

    function rhok = geometric(gamma,k)
        rhok = (RHO0 * gamma.^k)';
    end

end