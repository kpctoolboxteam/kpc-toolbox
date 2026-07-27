clear;

load("./data/24.hour.BuildServer.11-28-2007.07-55-PM.trace.csv.filtered.mat"); kpcsample = S(:,3);
kpcsample=kpcsample/mean(kpcsample);

options = kpcfit_ph_options(E,'MinNumStates',2,'MaxNumStates',8,'Runs',10);

[PH]=kpcfit_ph_auto(E,options); 

kpcfit_ph_summary


oldmethodsPH = PH;

for k=1:15
    E(k) = mean(kpcsample.^k);          
end

kpcsample=kpcsample/mean(kpcsample);

[F,xccdf] = ecdf(kpcsample);


xfiltered = xccdf(1:length(xccdf)/100:end);
xoldmethods = [(0:xfiltered(1)/10:xfiltered(1))';xfiltered];


figure;
loglog(xccdf,1-F); hold on;

labelsOldMethods{1} = 'trace';

indexFrom = 1;
indexTo = size(oldmethodsPH,1);

for i=indexFrom:indexTo
    oldmethods_cdf{i} = map_cdf(oldmethodsPH{i},xoldmethods);
    loglog(xoldmethods,1-oldmethods_cdf{i},"-."); hold on;

    method = oldmethodsPH{i,4};
    switch method
        case 'exact'
             labelsOldMethods{end+1} = strcat("Prony's - ",num2str(size(oldmethodsPH{i}{1},1)));
        case 'approx_moment_space'
             labelsOldMethods{end+1} = strcat("KPC - ",num2str(size(oldmethodsPH{i}{1},1)));
        case 'approx_param_space'
             labelsOldMethods{end+1} = strcat('HyperExp - ',num2str(size(oldmethodsPH{i}{1},1)));
        otherwise
             warning('Unexpected fitting method')
    end
end

legend(labelsOldMethods,Location='best')
title('CCDF of trace and fitted PHs')


figure;
plot(1:length(E),log(E),'o'); hold on;

for i=indexFrom:indexTo
    plot(1:length(E),log(map_moment(oldmethodsPH{i},1:length(E))),'x'); hold on;
end

legend(labelsOldMethods,Location='best')
title('Moments of trace and fitted PHs')


figure;
trunc = find(F<0.95);

plot(xccdf(trunc),F(trunc),'--',LineWidth=1.6); hold on;

for i=indexFrom:indexTo
    truncOldMethods = find(oldmethods_cdf{i}<0.95);
    plot(xoldmethods(truncOldMethods),oldmethods_cdf{i}(truncOldMethods),LineWidth=1.6); hold on;
end

xlabel('t - time') 
ylabel('F(t) - CDF') 
legend(labelsOldMethods,Location='best')
title('CDF of trace and fitted PHs')


save(trace+"_oldmethods.mat","oldmethodsPH","labelsOldMethods","kpcsample","oldmethods_cdf","xoldmethods")







