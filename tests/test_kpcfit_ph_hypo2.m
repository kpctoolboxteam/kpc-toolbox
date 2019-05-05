D0=[-2.389377046831244e+00     1.907112341413688e+00    
     1.141504016723399e-01    -4.480622953168757e+00]
D1=[ 0     4.822647054175562e-01
     4.349561509345400e+00     1.691104215101627e-02]
HYPO2={D0,D1};
E=map_moment(HYPO2,1:3);
options = kpcfit_ph_options(E);
fprintf('\n*******************************************************************************\n');
fprintf('Generating 3 moments from Hypo-exponential(2) (Generalized Erlang)\n')
MAP={D0,D1};
fprintf('Running kpc-ph fitting method\n')
PH = kpcfit_ph_auto(E, options);

for j = 1:length(PH)
    alpha = map_pie(PH{j})
    T = PH{j}{1}
end

fprintf('\n*******************************************************************************\n');
fprintf('Now setting E[X^3]=0 in Hypo-exponential(2) moments\n');
E(3)=0;
MAP={D0,D1};
fprintf('Running kpc-ph fitting method\n')
PH = kpcfit_ph_auto(E, options);

for j = 1:length(PH)
    alpha = map_pie(PH{j})
    T = PH{j}{1}
end
