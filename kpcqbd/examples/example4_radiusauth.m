setup_paths;

fprintf('=== KPC-QBD Fitting: RADIUS Authentication Trace ===\n\n');

datafile = '../data/simulation_probs_avgtable_radiusauth_1e10.mat';
if ~exist(datafile, 'file')
    error('Data file not found: %s\nPlease ensure the data file exists.', datafile);
end

load(datafile);
fprintf('Loaded simulation data from RADIUS Authentication trace\n');

rho_values = [0.25, 0.5, 0.75, 0.9];
rho_idx = 2;
rho = rho_values(rho_idx);

scv_arrival = 5;
nprobs = 20;

fprintf('\nSystem parameters:\n');
fprintf('  Utilization (rho): %.2f\n', rho);
fprintf('  Fitting %d queue length probabilities\n', nprobs);

simprobs = prob{rho_idx}(1:nprobs);
simprobs = simprobs(:);

ARV = hyperexp_create(1/rho, scv_arrival);

configurations = {
    struct('K', 2, 'J', 2, 'name', 'K=2, J=2 (small)'),
    struct('K', 3, 'J', 3, 'name', 'K=3, J=3 (medium)'),
    struct('K', 4, 'J', 3, 'name', 'K=4, J=3 (large)')
};

results = cell(length(configurations), 1);

fprintf('\n--- Fitting Multiple Configurations ---\n');

for c = 1:length(configurations)
    cfg = configurations{c};
    fprintf('\nConfiguration %d: %s\n', c, cfg.name);

    tic;
    [PH, score, eflag, x, PHj] = kpcqbd_fit(ARV, cfg.J, cfg.K, simprobs);
    elapsed = toc;

    fitted_probs = kpcqbd_solve(ARV, PHj, nprobs);
    fitted_probs = fitted_probs(:);
    rel_error = abs(fitted_probs - simprobs) ./ simprobs;

    results{c}.cfg = cfg;
    results{c}.score = score;
    results{c}.eflag = eflag;
    results{c}.time = elapsed;
    results{c}.fitted_probs = fitted_probs;
    results{c}.rel_error = rel_error;
    results{c}.PHj = PHj;
    results{c}.nstates = prod(cellfun(@(ph) length(ph{1}), PHj));

    fprintf('  Time: %.2f sec, Score: %.6f, States: %d\n', ...
        elapsed, score, results{c}.nstates);
    fprintf('  Mean rel. error: %.4f%%, Max rel. error: %.4f%%\n', ...
        mean(rel_error)*100, max(rel_error)*100);
end

fprintf('\n--- Configuration Comparison ---\n');
fprintf('%-20s  %8s  %8s  %8s  %8s\n', 'Configuration', 'Score', 'States', 'Mean Err', 'Time(s)');
fprintf('%-20s  %8s  %8s  %8s  %8s\n', '-------------', '-----', '------', '--------', '-------');
for c = 1:length(configurations)
    r = results{c};
    fprintf('%-20s  %8.4f  %8d  %7.2f%%  %8.2f\n', ...
        r.cfg.name, r.score, r.nstates, mean(r.rel_error)*100, r.time);
end

[~, best_idx] = min(cellfun(@(r) r.score, results));
best = results{best_idx};
fprintf('\nBest configuration: %s\n', best.cfg.name);

figure;

subplot(2,1,1);
semilogy(0:nprobs-1, simprobs, 'ko-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Simulation');
hold on;
colors = {[0.8 0.2 0.2], [0.2 0.6 0.2], [0.2 0.2 0.8]};
markers = {'^--', 's--', 'd--'};
for c = 1:length(configurations)
    semilogy(0:nprobs-1, results{c}.fitted_probs, markers{c}, ...
        'LineWidth', 1.5, 'MarkerSize', 6, 'Color', colors{c}, ...
        'DisplayName', configurations{c}.name);
end
hold off;
xlabel('Queue Length');
ylabel('Probability (log scale)');
title(sprintf('RADIUS Auth: Queue Length Distribution (\\rho = %.2f)', rho));
legend('Location', 'northeast');
grid on;

subplot(2,1,2);
hold on;
bar_width = 0.25;
for c = 1:length(configurations)
    x_pos = (0:nprobs-1) + (c-2)*bar_width;
    bar(x_pos, results{c}.rel_error * 100, bar_width, ...
        'FaceColor', colors{c}, 'DisplayName', configurations{c}.name);
end
hold off;
xlabel('Queue Length');
ylabel('Relative Error (%)');
title('Fitting Error Comparison');
legend('Location', 'northeast');
grid on;

fprintf('\n--- Sensitivity to Utilization ---\n');
fprintf('Using best configuration: %s\n\n', best.cfg.name);

figure;
colors_rho = {[0.2 0.4 0.8], [0.2 0.7 0.3], [0.9 0.5 0.1], [0.8 0.2 0.2]};

for r_idx = 1:length(rho_values)
    rho_curr = rho_values(r_idx);
    ARV_curr = hyperexp_create(1/rho_curr, scv_arrival);
    simprobs_curr = prob{r_idx}(1:nprobs);
    simprobs_curr = simprobs_curr(:);

    [~, ~, ~, ~, PHj_curr] = kpcqbd_fit(ARV_curr, best.cfg.J, best.cfg.K, simprobs_curr);
    fitted_curr = kpcqbd_solve(ARV_curr, PHj_curr, nprobs);
    fitted_curr = fitted_curr(:);

    rel_err = mean(abs(fitted_curr - simprobs_curr) ./ simprobs_curr) * 100;
    fprintf('  rho=%.2f: Mean relative error = %.2f%%\n', rho_curr, rel_err);

    subplot(2, 2, r_idx);
    semilogy(0:nprobs-1, simprobs_curr, 'ko-', 'LineWidth', 1.5, 'DisplayName', 'Sim');
    hold on;
    semilogy(0:nprobs-1, fitted_curr, 's--', 'Color', colors_rho{r_idx}, ...
        'LineWidth', 1.5, 'DisplayName', 'KPC-QBD');
    hold off;
    xlabel('Queue Length');
    ylabel('Probability');
    title(sprintf('\\rho = %.2f', rho_curr));
    legend('Location', 'northeast');
    grid on;
end
sgtitle('RADIUS Auth: KPC-QBD Fit Across Utilization Levels');

fprintf('\nExample 4 completed successfully.\n');
