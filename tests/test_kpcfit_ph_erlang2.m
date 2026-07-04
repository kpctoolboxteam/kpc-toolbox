E=map_moment(map_erlang(2,2),1:3);
options = kpcfit_ph_options(E);
fprintf('\n*******************************************************************************\n');
fprintf('Generating 3 moments from Erlang-2\n')
fprintf('Running kpc-ph fitting method\n')
PH = kpcfit_ph_auto(E,options);

for j = 1:size(PH,1)
    alpha = map_pie(PH{j,1})
    T = PH{j,1}{1}
end
