seed = 4;
%seed = 12;
%seed = 27;
rand('seed',seed)
load BCAUG89.mat
trace = kpcfit_init(S);
%%
MAP = kpcfit_auto(trace,'OnlyAC',1,'MaxRunsAC',1,'MaxIterAC',50,'AnimateAC',1)
%%
subplot(1,2,1)
plot(trace_acf(S,1:50),'k-'); hold on; plot(1:50,trace.ACFull(1:50),'r-');
plot(trace.ACLags,trace.AC,'r*');
plot(1:50,map_acf(MAP,1:50),'b-');
ylabel('Autocorrelation Function - \rho_k')
xlabel('Lag k')
xlim([1,50])

subplot(1,2,2)
loglog(trace_acf(S,1:max(trace.ACLags)),'k-'); hold on; plot(trace.ACLags,trace.AC,'r-');
plot(trace.ACLags,map_acf(MAP,trace.ACLags),'b-');
warning off
ylabel('Autocorrelation Function - \rho_k - [log]')
xlabel('Lag k - [log]')
warning on