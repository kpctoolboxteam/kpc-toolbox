%% Example 3: KPC-QBD Fitting with BuildServer Trace
% This example demonstrates KPC-QBD fitting using real simulation data
% from a build server workload trace.
%
% Prerequisites:
%   - KPC Toolbox must be on the MATLAB path
%   - Optimization Toolbox required
%   - Data file: ../data/simulation_probs_avgtable_buildserver_1e10.mat

%% Setup paths
setup_paths;

%% Load simulation data
fprintf('=== KPC-QBD Fitting: BuildServer Trace ===\n\n');

datafile = '../data/simulation_probs_avgtable_buildserver_1e10.mat';
if ~exist(datafile, 'file')
    error('Data file not found: %s\nPlease ensure the data file exists.', datafile);
end

load(datafile);  % Loads 'prob' cell array with simulation probabilities
fprintf('Loaded simulation data from BuildServer trace\n');
fprintf('  Simulation samples: 1e10\n');

%% Define parameters
% Available utilization levels: 0.25, 0.5, 0.75, 0.9
rho_values = [0.25, 0.5, 0.75, 0.9];
rho_idx = 1;                    % Use rho = 0.25
rho = rho_values(rho_idx);

scv_arrival = 5;                % SCV of arrival process (matches simulation)
nprobs = 20;                    % Number of probabilities to fit

% KPC parameters
K = 4;                          % Phase-type distribution size
J = 3;                          % Number of PH distributions

fprintf('\nSystem parameters:\n');
fprintf('  Utilization (rho): %.2f\n', rho);
fprintf('  Arrival SCV: %.1f\n', scv_arrival);
fprintf('  KPC configuration: K=%d, J=%d\n', K, J);

%% Extract target probabilities from simulation
simprobs = prob{rho_idx}(1:nprobs);
fprintf('\nSimulation probabilities (first 5):\n');
for i = 1:min(5, nprobs)
    fprintf('  P(N = %2d) = %.6f\n', i-1, simprobs(i));
end

%% Create arrival process
ARV = hyperexp_create(1/rho, scv_arrival);
fprintf('\nArrival process created:\n');
fprintf('  Mean: %.4f\n', 1/map_lambda(ARV));
fprintf('  SCV: %.4f\n', map_scv(ARV));

%% Run KPC-QBD fitting
fprintf('\nRunning KPC-QBD fitting...\n');
fprintf('  Optimizing %d parameters...\n', 2 + (J-1)*(K-1) + (J-1)*K);

tic;
[PH_result, score, eflag, x, PH_components, x0] = kpcqbd_fit(ARV, J, K, simprobs);
fitting_time = toc;

fprintf('  Fitting completed in %.2f seconds\n', fitting_time);
fprintf('  Optimization score: %.6f\n', score);
fprintf('  Exit flag: %d ', eflag);
if eflag > 0
    fprintf('(converged)\n');
else
    fprintf('(may not have converged)\n');
end

%% Compute fitted probabilities
fitted_probs = kpcqbd_solve(ARV, PH_components, nprobs);
fitted_probs = fitted_probs(:);  % Ensure column vector
simprobs = simprobs(:);          % Ensure column vector

%% Display fitted service process characteristics
fprintf('\nFitted service process (KPC composition of %d PHs):\n', J);
for j = 1:length(PH_components)
    ph = PH_components{j};
    fprintf('  PH%d: States=%d, Mean=%.4f, SCV=%.4f\n', ...
        j, length(ph{1}), 1/map_lambda(ph), map_scv(ph));
end

%% Compute fitting accuracy
rel_error = abs(fitted_probs - simprobs) ./ simprobs;
max_rel_error = max(rel_error);
mean_rel_error = mean(rel_error);

fprintf('\nFitting accuracy:\n');
fprintf('  Max relative error: %.4f%%\n', max_rel_error * 100);
fprintf('  Mean relative error: %.4f%%\n', mean_rel_error * 100);

%% Display probability comparison
fprintf('\nProbability comparison:\n');
fprintf('  State    Simulated   Fitted      Rel.Error\n');
fprintf('  -----    ---------   ------      ---------\n');
for i = 1:nprobs
    fprintf('  %3d     %.6f    %.6f     %.2f%%\n', ...
        i-1, simprobs(i), fitted_probs(i), rel_error(i)*100);
end

%% Plot results
figure;

% Subplot 1: Log-scale comparison (typical for queue length distributions)
subplot(2,1,1);
semilogy(0:nprobs-1, simprobs, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Simulation');
hold on;
semilogy(0:nprobs-1, fitted_probs, 'r^--', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'KPC-QBD Fit');
hold off;
xlabel('Queue Length');
ylabel('Probability (log scale)');
title(sprintf('BuildServer: Queue Length Distribution (\\rho = %.2f)', rho));
legend('Location', 'northeast');
grid on;

% Subplot 2: Relative error
subplot(2,1,2);
bar(0:nprobs-1, rel_error * 100, 'FaceColor', [0.4 0.6 0.8]);
xlabel('Queue Length');
ylabel('Relative Error (%)');
title('Fitting Error by Queue Length');
grid on;
if max(rel_error) > 0
    ylim([0, max(rel_error*100)*1.2]);
end

fprintf('\nExample 3 completed successfully.\n');
