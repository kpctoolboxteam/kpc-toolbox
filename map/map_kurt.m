function KURT=map_kurt(MAP)
% KURT=map_kurt(MAP) - Compute kurtosis
%
%  Input: 
%  MAP: a MAP in the form of {D0,D1}
%       
%  Output: 
%  KURT: kurtosis of the interarrrival times
%
%
m=zeros(1,4);
for i=1:4
    m(i)=map_moment(MAP,i);
end
KURT=(m(4)-4*m(3)*m(1)+6*m(2)*m(1)^2-3*m(1)^4)/map_var(MAP)^2;
