function MEAN=map_mean(MAP)
% MEAN=map_mean(MAP) - Compute mean inter-arrival time
%
%  Input:
%  MAP: a MAP in the form of {D0,D1}
%
%  Output:
%  MEAN: mean inter-arrival time
%

MEAN=1/map_lambda(MAP);
end