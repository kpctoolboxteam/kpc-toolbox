seed = 1;
rand('seed',seed)
load BCAUG89.mat
trace = kpcfit_init(S);
MAP = kpcfit_auto(trace,'OnlyAC',1,'MaxIterAC',50);
BCAUG_ORD2 = MAP;

%%
MAP = kpcfit_auto(trace,'MaxRunsAC',10,'MaxIterBC',5,'MaxRunsBC',1)

% %%
% subplot(1,2,1)
% plot(trace_acf(S,1:50),'k-'); hold on; plot(1:50,trace.ACFull(1:50),'r-');
% plot(trace.ACLags,trace.AC,'r*');
% plot(1:50,map_acf(MAP,1:50),'b-');
% ylabel('Autocorrelation Function - \rho_k')
% xlabel('Lag k')
% xlim([1,50])
% 
% subplot(1,2,2)
% loglog(trace_acf(S,1:max(trace.ACLags)),'k-'); hold on; plot(trace.ACLags,trace.AC,'r-');
% plot(trace.ACLags,map_acf(MAP,trace.ACLags),'b-');
% ylabel('Autocorrelation Function - \rho_k - [log]')
% xlabel('Lag k - [log]')


%%
figure

subplot(1,2,2)
BCest=[];
surf(log10(trace.BCGridLags)',log10(trace.BCGridLags)',reshape((trace.BC),5,5)); hold all;
for i=1:size(trace.BCLags,1)
    BCest(end+1) = map_joint(MAP,trace.BCLags(i,:),[1,1,1]);
end
surf(log10(trace.BCGridLags)',log10(trace.BCGridLags)',reshape((BCest),5,5)); hold all;
title('third-order fit')

subplot(1,2,1)
BCest2=[];
surf(log10(trace.BCGridLags)',log10(trace.BCGridLags)',reshape((trace.BC),5,5)); hold all;
for i=1:size(trace.BCLags,1)
    BCest2(end+1) = map_joint(BCAUG_ORD2,trace.BCLags(i,:),[1,1,1]);
end
surf(log10(trace.BCGridLags)',log10(trace.BCGridLags)',reshape((BCest2),5,5)); hold all;
title('second-order fit')
