function PIQ=map_prob(MAP)
% PIQ=map_prob(MAP) - Equilibrium distribution of the underlying
% continuous-time process
%
%  Input:
%  MAP: a MAP in the form of {D0,D1}
%
%  Output:
%  PIQ: equilibrium distribution of the continuous-time Markov chain
%  Q=D0+D1
%

PIQ=ctmc_solve(map_infgen(MAP)); 
end
