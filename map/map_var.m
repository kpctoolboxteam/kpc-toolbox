function VAR=map_var(MAP)
% VAR=map_var(MAP) - Compute variance of interarrival times
%
%  Input:
%  MAP: a MAP in the form of {D0,D1}
%
%  Output:
%  VAR: variance of inter-arrival times
%

VAR=map_moment(MAP,2)-map_mean(MAP)^2;
end