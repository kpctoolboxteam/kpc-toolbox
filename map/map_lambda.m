function LAMBDA=map_lambda(MAP)
% LAMBDA=map_lambda(MAP) - Compute mean arrival rate
%
%  Input:
%  MAP: a MAP in the form of {D0,D1}
%
%  Output:
%  LAMBDA: mean arrival rate (1/mean inter-arrival time)
%
% MAP Queueing Networks Toolbox
% Version 1.0 	 15-Apr-2008
p=map_prob(MAP);
LAMBDA=p*MAP{2}*ones(size(MAP{2},2),1);
end