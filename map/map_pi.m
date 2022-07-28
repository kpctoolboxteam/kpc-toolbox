function p=map_pi(MAP)
% MAP_PI computes the stationary probability distribution of the phase in 
% a MAP, at arrival moments
% MAP:  Markovian arrival process to evaluate
%
% Copyright (c) 2012-2018, Imperial College London 
% All rights reserved.

% p=dtmc_solve(map_embedded(MAP)); % this does not work with
% E1=map_erlang(1,10), the first row becomes perfectly 0 

p = map_pie(MAP);
end
