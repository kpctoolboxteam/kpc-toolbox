function MAP=map_exponential(MEAN)
% MAP=map_exponential(MEAN) - Fit a Poisson process as a MAP
%
%  Input: 
%  MEAN: mean inter-arrival time of the process
%       
%  Output: 
%  MAP: a MAP in the form of {D0,D1}
%
%  Examples:
%  - MAP=map_exponential(2) return a Poisson process with rate lambda=0.5
%

if iscell(MEAN)
    MAP=map_exponential(map_mean(MEAN));
    return
end
mu=1/MEAN;
MAP={[-mu],[mu]};
end