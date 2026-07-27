function MAPOUT=map_renewal(MAPIN)
% MAPOUT=map_renewal(MAPIN) - Remove all correlations
%
%  Input:
%  MAPIN: a MAP in the form of {D0,D1}
%
%  Output:
%  MAPOUT: renewal MAP process with same cdf as MAPIN, but no
%  correlations between inter-arrival times
%

MAPOUT=MAPIN;
MAPOUT{2}=MAPIN{2}*ones((length(MAPIN{2})),1)*map_pie(MAPIN);
end