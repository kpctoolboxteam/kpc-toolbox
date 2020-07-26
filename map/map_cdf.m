function CDFVALS=map_cdf(MAP,POINTS)
% CDFVALS=map_cdf(MAP,POINTS) - Compute cumulative distribution function of
% interarrival times
%
%  Input: 
%  MAP: a MAP in the form of {D0,D1}
%  POINTS: a set of interarrival times
%       
%  Output: 
%  CDFVALS: values of the cumulative distribution at the specified points
%  returned in the same order of the POINTS vector
%
%  Examples:
%  - map_cdf(MAP,1) returns Pr[T<=1], being T an interarrival time
%  - map_cdf(MAP,[1,5]) returns [Pr[T<=1], Pr[T<=5]]

CDFVALS=0*POINTS(:)';
pie = map_pie(MAP);
e1=ones(length(MAP{2}),1);
for t=1:length(POINTS)
    CDFVALS(t)=1-pie*expm(MAP{1}*POINTS(t))*e1;
end
end