%% plot distance function
dist_exact = [];
ord_exact = [];
mom_exact = [];
best_exact = Inf;
best_exact_pos = 0;
dist_apxms = [];
ord_apxms = [];
mom_apxms = [];
best_apxms = Inf;
best_apxms_pos = 0;
dist_apxps = [];
ord_apxps = [];
mom_apxps = [];
best_apxps = Inf;
best_apxps_pos = 0;
for i=1:size(PH,1)
    if strcmpi(PH{i,4},'exact')
        ord_exact(end+1) = length(PH{i,1}{1});
        dist_exact(end+1) = PH{i,2};
        mom_exact(end+1,1:length(E)) = map_moment(PH{i,1},1:length(E));
        if best_exact_pos == 0 || best_exact > dist_exact(end)
            best_exact = dist_exact(end);
            best_exact_pos = i;
            best_exact_relpos = length(ord_exact);
        end
    elseif strcmpi(PH{i,4},'approx_param_space')
        ord_apxps(end+1) = length(PH{i,1}{1});
        dist_apxps(end+1) = PH{i,2};
        mom_apxps(end+1,1:length(E)) = map_moment(PH{i,1},1:length(E));        
        if best_apxps_pos == 0 || best_apxps > dist_apxps(end)
            best_apxps = dist_apxps(end);
            best_apxps_pos = i;
            best_apxps_relpos = length(ord_apxps);
        end
    elseif strcmpi(PH{i,4},'approx_moment_space')
        ord_apxms(end+1) = length(PH{i,1}{1});
        dist_apxms(end+1) = PH{i,2};
        mom_apxms(end+1,1:length(E)) = map_moment(PH{i,1},1:length(E));        
        if best_apxms_pos == 0 || best_apxms > dist_apxms(end)
            best_apxms = dist_apxms(end);
            best_apxms_pos = i;
            best_apxms_relpos = length(ord_apxms);
        end
    end
end
%% 
% loglog(ord_exact,dist_exact,'x-');
% hold all;
% loglog(ord_apxms,dist_apxms,'o-');
% loglog(ord_apxps,dist_apxps,'*-');
% 
% ylabel('moment set distance')
% xlabel('number of states (PH order)')
% xlim([min([ord_exact,ord_apxms,ord_apxps]/2),max([ord_exact,ord_apxms,ord_apxps]*2)])
% legend(sprintf('best exact mom. match (idx=%d)',best_exact_pos), sprintf('best approx - mom. space (idx=%d)',best_apxps_pos), sprintf('best approx - param. space (idx=%d)',best_apxms_pos), 'Location', 'Best')
% title(trace_name)
%% moment values
fprintf('Log10\tTRACE        ',j)
for j=1:length(ord_exact)
    fprintf('\tEXMM (idx=%d)   ',j)
end
for j=1:length(ord_apxms)
    fprintf('\tAPXMS (idx=%d)  ',length(ord_exact)+j)
end
for j=1:length(ord_apxps)
    fprintf('\tAPXPS (idx=%d)  ',length(ord_exact)+length(ord_apxms)+j)
end
for k=1:length(E)
    if k==1
        fprintf('\nE[X]')
    else
        fprintf('\nE[X^%d]',k)
    end
    fprintf('\t%13d',log10(E(k)))
    for j=1:length(ord_exact)
        fprintf('\t%13d',log10(mom_exact(j,k)))
    end
    for j=1:length(ord_apxms)
        fprintf('\t%13d',log10(mom_apxms(j,k)+eps))
    end
    for j=1:length(ord_apxps)
        fprintf('\t%13d',log10(mom_apxps(j,k)+eps))
    end
end
fprintf('\n')
%% best moment plot
figure
plot(1:length(E),log(E),'k-')
hold all
if exist('best_exact_relpos')
plot(1:length(E),log(mom_exact(best_exact_relpos,1:length(E))),'x')
end
if exist('best_apxms_relpos')
plot(1:length(E),log(mom_apxms(best_apxms_relpos,1:length(E))),'o')
end
if exist('best_apxps_relpos')
plot(1:length(E),log(mom_apxps(best_apxps_relpos,1:length(E))),'d')
end
ylabel('log E[X^k]')
xlabel('moment order k')
legend('trace', sprintf('best exact mom. match (idx=%d)',best_exact_pos), sprintf('best approx - mom. space (idx=%d)',best_apxms_pos), sprintf('best approx - param. space (idx=%d)',best_apxps_pos), 'Location', 'Best')
%title(trace_name)

