function [GAMMA, RHO0, RESIDUALS] = map_gamma(MAP, limit)
% Estimates the auto-correlation decay rate of a MAP.
% For MAPs of order higher than 2, performs an approximation of the ACF
% curve using non-linear least squares fitting.
% Input:
% - MAP: the MAP
% - limit: maximum lag considered (optional, default = 1000)
% Output:
% - GAMMA: autocorrelation decay rate

if nargin < 2
    limit = 1000;
    lag = 1:(limit/10):limit;    
end

if length(limit) ~= 1
	error('Invalid parameter'); 
end

n = size(MAP{1}, 1);

if n == 1
    % poisson process: no correlation
    GAMMA = 0;
elseif n == 2
    % second-order MAP: geometric ACF
    if abs(map_acf(MAP,1)) < 1e-8
        % phase-type
        GAMMA = 0;
    else
        GAMMA = map_acf(MAP,2) / map_acf(MAP,1);
    end
else
    % higher-order MAP: non-geometric
    
    M1 = map_mean(MAP);
    M2 = map_moment(MAP, 2);
    VAR = M2-M1^2;
    SCV = VAR/M1^2;
    RHO0 = 1/2 * (1 - 1/SCV);
    
    rho = map_acf(MAP, lag)';
    
    %problem.Variables = 1;
    %problem.LB = -1;
    %problem.UB = 0.999;
    %problem.ObjFunction = @(x) sum((geometric(x,lag)-rho).^2);
    %GAMMA = PSwarm(problem);

    opt = statset('nlinfit');
    opt.MaxIter = 1e5;
    opt.Display = 'off';
    opt.RobustWgtFun = 'fair';
    try
        [GAMMA,RESIDUALS] = nlinfit(lag, rho, @geometric, 0.99, opt);
    catch me
        warning('Non linear regression for ACF decay rate failed, trying lsqcurvefit'); 
        GAMMA = lsqcurvefit(@geometric, 0.99, lag, rho, -1, 1);
    end
end

    function rhok = geometric(gamma,k)
        rhok = (RHO0 * gamma.^k)';
    end

end