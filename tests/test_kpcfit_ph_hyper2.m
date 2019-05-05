H20{1}=[
    -12.000975230358000   0.000975230358000
    0.000080862578000  -0.088000599518000
    ];
H20{2}=[
    12.000000000000000                   0
    0   0.087919736940000
    ];
H20=map_renewal(H20);

E = map_moment(H20,1:3);
options = kpcfit_ph_options(E);
fprintf('\n*******************************************************************************\n');
fprintf('Generating 3 moments from Hyper-Exponential(3)\n')
MAP={D0,D1};
fprintf('Running kpc-ph fitting method\n')
PH = kpcfit_ph_auto(E,options);

for j = 1:length(PH)
    alpha = map_pie(PH{j})
    T = PH{j}{1}
end
