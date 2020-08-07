function MOMENTS=map_moment(MAP,ORDERS)
% MOMENTS=map_moment(MAP,ORDERS) -  Compute (power) moments of interarrival
% times of the specified order
%
%  Input:
%  MAP: a MAP in the form of {D0,D1}
%  ORDERS: set of moment orders (1=>E[X], 2=>E[X^2], ...)
%
%  Output:
%  MOMENTS: moments returned in the same order of ORDERS
%
%  Examples:
%  - map_moment(MAP,1:2) return the first two power moments E[X], E[X^2]
%  - map_moment(MAP,2) return the second power moment E[X^2]
%
if MAP{1}==0
    MOMENTS = 0*ORDERS;
    return
end

D0=MAP{1};
D1=MAP{2};
e=ones(length(D0),1);
%pi=map_prob(MAP);
%x=pi*D1;
%y=e/(pi*D1*e);
x=map_pie(MAP);
y=e;
for t=1:length(ORDERS)
    i=ORDERS(t);
    if isnan(D0)
        MOMENTS(t)=NaN;
    else
        A=(inv(-D0)^i);
        MOMENTS(t)=factorial(i)*x*A*y;
    end
end
end
