%% Example 2: KPC-QBD Fitting with Synthetic Target
% This example demonstrates the KPC-QBD fitting algorithm. We:
%   1. Generate target queue length probabilities from a known distribution
%   2. Fit a KPC-QBD model to match these probabilities
%   3. Compare fitted vs target distributions
%
% Prerequisites:
%   - KPC Toolbox must be on the MATLAB path
%   - Optimization Toolbox required

%% Setup paths
setup_paths;

%% Define system parameters
rho = 0.5;              % System utilization
scv_arrival = 5;        % SCV for arrival process
nprobs = 20;            % Number of probabilities to fit

%% Fitting parameters
K = 3;                  % Phase-type distribution size
J = 2;                  % Number of PH distributions in KPC

fprintf('=== KPC-QBD Fitting Example ===\n');
fprintf('Parameters: rho=%.2f, K=%d, J=%d\n\n', rho, K, J);

%% Create arrival process (same for target and fitted model)
ARV = hyperexp_create(1/rho, scv_arrival);

%% Generate target probabilities using a known service distribution
fprintf('Generating target queue length probabilities...\n');

% Create a "ground truth" service process
% Using manually constructed phase-type distributions

% Target PH1: Erlang-2 like (SCV < 1)
mu1 = 2;
target_PH1 = {[-mu1, mu1; 0, -mu1], [0, 0; mu1, 0]};

% Target PH2: Hyperexponential (SCV > 1)
p2 = 0.4;
lambda1 = 3; lambda2 = 0.6;
target_PH2 = {diag([-lambda1, -lambda2]), [lambda1; lambda2] * [p2, 1-p2]};

target_H = {target_PH1, target_PH2};
target_probs = kpcqbd_solve(ARV, target_H, nprobs);

fprintf('  Target distribution generated (%d probabilities)\n', nprobs);
fprintf('  Target P(empty) = %.4f\n', target_probs(1));

%% Run KPC-QBD fitting
fprintf('\nRunning KPC-QBD fitting optimization...\n');
fprintf('  This may take a few minutes...\n');

tic;
[PH_fitted, score, eflag, x, PH_components] = kpcqbd_fit(ARV, J, K, target_probs);
fitting_time = toc;

fprintf('  Fitting completed in %.2f seconds\n', fitting_time);
fprintf('  Optimization score (relative error): %.6f\n', score);
fprintf('  Exit flag: %d\n', eflag);

%% Compute fitted probabilities
fitted_probs = kpcqbd_solve(ARV, PH_components, nprobs);

%% Display fitted PH characteristics
fprintf('\nFitted PH components:\n');
for j = 1:length(PH_components)
    ph = PH_components{j};
    fprintf('  PH%d: Mean=%.4f, SCV=%.4f, States=%d\n', ...
        j, 1/map_lambda(ph), map_scv(ph), length(ph{1}));
end

%% Compute error metrics
abs_error = abs(fitted_probs - target_probs);
rel_error = abs_error ./ target_probs;
max_abs_error = max(abs_error);
max_rel_error = max(rel_error);
mean_rel_error = mean(rel_error);

fprintf('\nFitting accuracy:\n');
fprintf('  Max absolute error: %.6f\n', max_abs_error);
fprintf('  Max relative error: %.4f%%\n', max_rel_error * 100);
fprintf('  Mean relative error: %.4f%%\n', mean_rel_error * 100);

%% Display probability comparison
fprintf('\nProbability comparison (first 10 states):\n');
fprintf('  State    Target      Fitted      Rel.Error\n');
fprintf('  -----    ------      ------      ---------\n');
for i = 1:min(10, nprobs)
    fprintf('  %3d     %.6f    %.6f     %.2f%%\n', ...
        i-1, target_probs(i), fitted_probs(i), rel_error(i)*100);
end

%% Plot comparison
figure;

% Subplot 1: Probability distributions
subplot(2,1,1);
bar_width = 0.35;
x_pos = 0:nprobs-1;
bar(x_pos - bar_width/2, target_probs, bar_width, 'FaceColor', [0.3 0.5 0.7], 'DisplayName', 'Target');
hold on;
bar(x_pos + bar_width/2, fitted_probs, bar_width, 'FaceColor', [0.8 0.4 0.3], 'DisplayName', 'Fitted');
hold off;
xlabel('Queue Length');
ylabel('Probability');
title(sprintf('Queue Length Distribution Comparison (\\rho = %.2f)', rho));
legend('Location', 'northeast');
grid on;

% Subplot 2: Relative error
subplot(2,1,2);
bar(x_pos, rel_error * 100, 'FaceColor', [0.5 0.7 0.5]);
xlabel('Queue Length');
ylabel('Relative Error (%)');
title('Fitting Error by Queue Length');
grid on;

fprintf('\nExample 2 completed successfully.\n');
