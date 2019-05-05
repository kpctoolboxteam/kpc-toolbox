function PROB=dtmc_solve(P,options)
% PROB=dtmc_solve(P) - Equilibrium distribution of a discrete-time Markov
% chain
%
%  Input:
%  P: stochastic transition matrix of the discrete-time Markov chain
%
%  Output:
%  PROB: equilibrium distribution of the discrete-time Markov chain
%
%  Examples:
%  - dtmc_solve([0.5,0.5;0.2,0.8])
%

if nargin<2
    PROB=ctmc_solve(P-eye(size(P)));
else
    PROB=ctmc_solve(P-eye(size(P)),options);
end
end