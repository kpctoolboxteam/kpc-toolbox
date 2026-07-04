function ARV = hyperexp_create(mean_val, scv)
%HYPEREXP_CREATE Create hyperexponential MAP representation
%   ARV = HYPEREXP_CREATE(MEAN, SCV) creates a 2-phase hyperexponential
%   distribution with given mean and squared coefficient of variation.
%
%   This is a standalone function that doesn't require Java dependencies.
%
%   Example:
%       ARV = hyperexp_create(2, 5);  % Mean=2, SCV=5

if scv <= 1
    error('SCV must be > 1 for hyperexponential distribution');
end

% Balanced hyperexponential parameters
p = 0.5 * (1 + sqrt((scv - 1) / (scv + 1)));
mu1 = 2 * p / mean_val;
mu2 = 2 * (1 - p) / mean_val;

% MAP representation {D0, D1}
D0 = diag([-mu1, -mu2]);
D1 = [mu1; mu2] * [p, 1-p];

ARV = {D0, D1};
end
