function VAR=map_var(MAP)
% VAR=map_var(MAP) - Compute variance of interarrival times
%
%  Input:
%  MAP: a MAP in the form of {D0,D1}
%
%  Output:
%  VAR: variance of inter-arrival times
%
% MAP Queueing Networks Toolbox
% Version 1.0 	 15-Apr-2008
VAR=map_moment(MAP,2)-map_mean(MAP)^2;
end