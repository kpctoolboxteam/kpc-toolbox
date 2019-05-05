function P=map_embedded(MAP)
% P=map_embedded(MAP) - Embedded discrete-time process
%
%  Input: 
%  MAP: a MAP in the form of {D0,D1}
%       
%  Output: 
%  P: transition matrix of the embedded process. If the MAP is feasible,
%  then P must be an irreducible stochastic matrix.
%
%  Examples:
%  - P=map_embedded(MAP) gives the probabilities P(i,j) that the MAP
%  restarts in phase j if the last absorption occurred in phase i
%

P=inv(-MAP{1})*MAP{2};
end