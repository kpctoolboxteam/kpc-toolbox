function MAPOUT=map_scale(MAPIN,NEWMEAN)
% MAPOUT=map_scale(MAPIN,NEWMEAN) - Rescale mean inter-arrival time
%
%  Input:
%  MAPIN: a MAP in the form of {D0,D1}
%  NEWMEAN: new mean inter-arrival time
%
%  Output:
%  MAPOUT: MAP with same normalized moments and correlations of MAPIN except
%  for the mean inter-arrival time that is set to NEWMEAN
%
%  Examples:
%  - map_mean(map_scale(map_exponential(1),2)) has mean 2
%

ratio=map_mean(MAPIN)/NEWMEAN;
MAPIN{1}=MAPIN{1}*ratio;
MAPIN{2}=MAPIN{2}*ratio;
MAPOUT=map_normalize(MAPIN);
end