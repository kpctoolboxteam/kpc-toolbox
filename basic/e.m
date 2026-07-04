function n=e(t)
% n = e(t) - column vector of ones.
%   e()      returns [1;1]
%   e(k)     returns ones(k,1)
%   e(MAP)   returns ones(order,1) for a MAP/cell argument
if nargin<1
    n=[1;1];
elseif iscell(t)
    n=ones(length(t{1}),1);
else
    n=ones(t,1);
end
end
