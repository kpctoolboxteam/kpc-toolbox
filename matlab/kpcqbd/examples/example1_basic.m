setup_paths;

rho = 0.5;
scv_arrival = 5;
nprobs = 30;

ARV = hyperexp_create(1/rho, scv_arrival);

fprintf('Arrival process:\n');
fprintf('  Mean inter-arrival time: %.4f\n', 1/map_lambda(ARV));
fprintf('  SCV: %.4f\n', map_scv(ARV));


mu1 = 2;
D0_ph1 = [-mu1, mu1; 0, -mu1];
D1_ph1 = [0, 0; mu1, 0];
PH1 = {D0_ph1, D1_ph1};

probs = [0.3, 0.5, 0.2];
rates = [2, 1, 0.5];
D0 = diag(-rates);
D1 = rates' * probs;
PH2 = {D0, D1};

fprintf('\nService process components:\n');
fprintf('  PH1 - Mean: %.4f, SCV: %.4f\n', 1/map_lambda(PH1), map_scv(PH1));
fprintf('  PH2 - Mean: %.4f, SCV: %.4f\n', 1/map_lambda(PH2), map_scv(PH2));

H = {PH1, PH2};

fprintf('\nSolving QBD...\n');
pa = kpcqbd_solve(ARV, H, nprobs);

fprintf('\nQueue length probabilities (first 10 states):\n');
for i = 1:min(10, nprobs)
    fprintf('  P(N = %2d) = %.6f\n', i-1, pa(i));
end

fprintf('\nSummary statistics:\n');
fprintf('  P(empty) = %.4f\n', pa(1));
fprintf('  Expected queue length = %.4f\n', (0:nprobs-1) * pa);
fprintf('  Sum of probabilities (first %d) = %.6f\n', nprobs, sum(pa));

figure;
bar(0:nprobs-1, pa, 'FaceColor', [0.3 0.5 0.7]);
xlabel('Queue Length');
ylabel('Probability');
title(sprintf('Queue Length Distribution (\\rho = %.2f)', rho));
grid on;

fprintf('\nExample 1 completed successfully.\n');
