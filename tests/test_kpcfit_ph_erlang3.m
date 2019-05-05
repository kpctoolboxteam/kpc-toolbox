E=map_moment(map_erlang(3,3),1:3);
options = kpcfit_ph_options(E);
fprintf('\n*******************************************************************************\n');
fprintf('Generating 3 moments from Erlang-3\n')
fprintf('Running kpc-ph fitting method\n')
PH = kpcfit_ph_auto(E,options);

for j = 1:length(PH)
    alpha = map_pie(PH{j})
    T = PH{j}{1}
end

%%
E=map_moment(map_erlang(3,3),1:3);
options = kpcfit_ph_options(E);
fprintf('\n*******************************************************************************\n');
fprintf('Generating 5 moments from Erlang-3\n')
fprintf('Running kpc-ph fitting method\n')
PH = kpcfit_ph_auto(E,options);

for j = 1:length(PH)
    alpha = map_pie(PH{j})
    T = PH{j}{1}
end
