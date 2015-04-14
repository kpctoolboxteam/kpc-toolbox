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
% MAP Queueing Networks Toolbox
% Version 1.0 	 15-Apr-2008
%
A=map_prob(MAP)*MAP{2};
PIE=A/(A*e(length(MAP{2})));
end
