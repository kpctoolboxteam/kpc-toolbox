%S=S/mean(S); 
for k=1:11 E(k) = mean(S.^k); end
options = kpcfit_ph_options(E,'MinNumStates',min(states),'MaxNumStates',max(states),'Runs',maxruns);
PH = kpcfit_ph_auto(E,options);
