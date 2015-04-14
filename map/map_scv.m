function SCV=map_scv(MAP)
% SCV=map_scv(MAP) - Compute squared coefficient of variation of
% inter-arrival times
%
%  Input:
%  MAP: a MAP in the form of {D0,D1}
%
%  Output:
%  SCV: squared coefficient of variation of inter-arrival times
%
% MAP Queueing Networks Toolbox
% Version 1.0 	 15-Apr-2008
SCV=map_var(MAP)/map_mean(MAP)^2;
end
