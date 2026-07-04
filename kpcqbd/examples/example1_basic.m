%% Example 1: Basic KPC-QBD Setup and Solution
% This example demonstrates the basic usage of kpcqbd_solve to compute
% steady-state queue length probabilities for a simple queueing system.
%
% Prerequisites:
%   - KPC Toolbox must be on the MATLAB path
%
% The example creates a simple G/G/1 queue with:
%   - HyperExponential arrival process
%   - Service process composed of two phase-type distributions

%% Setup paths
setup_paths;

%% Define system parameters
rho = 0.5;          % System utilization
scv_arrival = 5;    % Squared coefficient of variation for arrivals
nprobs = 30;        % Number of queue length probabilities to compute

%% Create arrival process
% HyperExponential with mean 1/rho and SCV = scv_arrival
ARV = hyperexp_create(1/rho, scv_arrival);

fprintf('Arrival process:\n');
fprintf('  Mean inter-arrival time: %.4f\n', 1/map_lambda(ARV));
fprintf('  SCV: %.4f\n', map_scv(ARV));

%% Create service process using KPC composition
% We create a service process as a Kronecker product composition of
% two simpler phase-type distributions

% First PH: Erlang-2 (low variability, SCV = 0.5)
% Two sequential exponential phases with rate mu each
mu1 = 2;  % rate for each phase
D0_ph1 = [-mu1, mu1; 0, -mu1];
D1_ph1 = [0, 0; mu1, 0];
PH1 = {D0_ph1, D1_ph1};

% Second PH: Categorical distribution (hyperexponential-like)
% 3-state categorical with different rates
probs = [0.3, 0.5, 0.2];
rates = [2, 1, 0.5];
D0 = diag(-rates);
D1 = rates' * probs;
PH2 = {D0, D1};

fprintf('\nService process components:\n');
fprintf('  PH1 - Mean: %.4f, SCV: %.4f\n', 1/map_lambda(PH1), map_scv(PH1));
fprintf('  PH2 - Mean: %.4f, SCV: %.4f\n', 1/map_lambda(PH2), map_scv(PH2));

% Combine using KPC
H = {PH1, PH2};

%% Solve the QBD
fprintf('\nSolving QBD...\n');
pa = kpcqbd_solve(ARV, H, nprobs);

%% Display results
fprintf('\nQueue length probabilities (first 10 states):\n');
for i = 1:min(10, nprobs)
    fprintf('  P(N = %2d) = %.6f\n', i-1, pa(i));
end

fprintf('\nSummary statistics:\n');
fprintf('  P(empty) = %.4f\n', pa(1));
fprintf('  Expected queue length = %.4f\n', (0:nprobs-1) * pa);
fprintf('  Sum of probabilities (first %d) = %.6f\n', nprobs, sum(pa));

%% Plot queue length distribution
figure;
bar(0:nprobs-1, pa, 'FaceColor', [0.3 0.5 0.7]);
xlabel('Queue Length');
ylabel('Probability');
title(sprintf('Queue Length Distribution (\\rho = %.2f)', rho));
grid on;

fprintf('\nExample 1 completed successfully.\n');
