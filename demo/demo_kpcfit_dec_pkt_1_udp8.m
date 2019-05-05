seed = 3
rand('seed',seed);
load DEC-PKT-1-UDP.mat
trace = kpcfit_init(S,'Smooth',2);
%%
[MAP,fac,fbc,kpcMAPs] = kpcfit_auto(trace,'OnlyAC',1,'MaxRunsAC',1,'MaxIterAC',50,'AnimateAC',1,'NumStates',8);
%%
subplot(1,2,1)
plot(trace_acf(S,1:100),'k-'); hold on; plot(1:100,trace.ACFull(1:100),'r-');
plot(trace.ACLags,trace.AC,'r*');
plot(1:100,map_acf(MAP,1:100),'b-');
ylabel('Autocorrelation Function - \rho_k')
xlabel('Lag k')
xlim([1,100])

subplot(1,2,2)
loglog(trace_acf(S,1:max(trace.ACLags)),'k-'); hold on; plot(trace.ACLags,trace.AC,'r-');
plot(trace.ACLags,map_acf(MAP,trace.ACLags),'b-');
ylabel('Autocorrelation Function - \rho_k - [log]')
xlabel('Lag k - [log]')
