function MEAN=map_mean(MAP)
% MEAN=map_mean(MAP) - Compute mean inter-arrival time
%
%  Input:
%  MAP: a MAP in the form of {D0,D1}
%
%  Output:
%  MEAN: mean inter-arrival time
%
% MAP Queueing Networks Toolbox
% Version 1.0 	 15-Apr-2008
MEAN=map_moment(MAP,1);
end