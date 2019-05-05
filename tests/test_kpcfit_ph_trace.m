function PH=test_kpcfit_ph_trace(trace_name,seed,states,maxruns)
fprintf(sprintf('Loading %s trace and generating moments\n\n', trace_name))
load(sprintf('%s.mat',trace_name))
if ~isempty(seed)
    rand('seed',seed);
end
%%
test_kpcfit_ph_run
%%
test_kpcfit_ph_summary
%%
figure
cdf_idx=round(linspace(1,length(X),min([1000,length(X)])));
loglog(X(cdf_idx),1-F(cdf_idx),'k-');
hold all
if best_exact_pos > 0
    loglog(X(cdf_idx),1-map_cdf(PH{best_exact_pos,1},X(cdf_idx)'),'-');
end
if best_apxms_pos > 0
    loglog(X(cdf_idx),1-map_cdf(PH{best_apxms_pos,1},X(cdf_idx)'),'-');
end
if best_apxps_pos > 0
    loglog(X(cdf_idx),1-map_cdf(PH{best_apxps_pos,1},X(cdf_idx)'),'-');
end
xlabel('Inter-arrival time - t');
ylabel('CCDF - Pr( X > t)');
if best_exact_pos > 0 && best_apxms_pos > 0 && best_apxps_pos > 0
    legend('trace',sprintf('best exact mom. match (idx=%d)',best_exact_pos), sprintf('best approx - mom. space (idx=%d)',best_apxms_pos), sprintf('best approx - param. space (idx=%d)',best_apxps_pos), 'Location', 'Best')
elseif best_apxms_pos > 0 && best_apxps_pos > 0
    legend('trace', sprintf('best approx - mom. space (idx=%d)',best_apxms_pos), sprintf('best approx - param. space (idx=%d)',best_apxps_pos), 'Location', 'Best')
end
title(trace_name)
end