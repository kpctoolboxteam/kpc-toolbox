function pos=minpos(v,n)
% pos = MINPOS(v)
% Position of maxima in vector v
% pos = MAXPOS(v,n)
% Position of n largest elements in vector v
%
% Copyright (c) 2012-2020, Imperial College London
% All rights reserved.
if nargin<2, n=1; end

if n==1
    [~,pos] = min(v);
else
    [~,pos] = sort(v, 'ascend');
    pos = pos(1:n);
end
    
end