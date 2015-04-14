function SKEWNESS=map_skew(MAP)
% SKEWNESS=map_skew(MAP) - Compute skewness of inter-arrival times
%
%  Input:
%  MAP: a MAP in the form of {D0,D1}
%
%  Output:
%  SKEWNESS: skewness of inter-arrival times
%

for i=1:3
    m(i)=map_moment(MAP,i);
end
M3=m(3)-3*m(2)*m(1)+2*m(1)^3;
SKEWNESS=M3/(sqrt(map_scv(MAP))*m(1))^3;
end
