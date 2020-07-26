function PIE=map_pie(MAP)
% PIE=map_pie(MAP) - Equilibrium distribution of the embedded
% discrete-time process
%
%  Input:
%  MAP: a MAP in the form of {D0,D1}
%
%  Output:
%  PIE: equilibrium distribution of the discrete-time Markov chain embedded
%  at departure instants P=map_embedded(MAP)=((-D0)^-1)*D1
%

A=map_prob(MAP)*MAP{2};
PIE=A/(A*ones(length(MAP{2}),1));
end
