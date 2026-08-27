function pos=minpos(v,n)
% pos = MINPOS(v)
% Position of the minimum of vector v (1-based, as in MIN)
% pos = MINPOS(v,n)
% Positions of the n smallest elements of vector v, ascending (1-based)
%
% NOTE: the Python kpctoolbox.minpos returns 0-based positions, following
% Python indexing. The values are the same, the offset is not.
%
% Copyright (c) 2012-2022, Imperial College London
% All rights reserved.
if nargin<2, n=1; end

if n==1
    [~,pos] = min(v);
else
    [~,pos] = sort(v, 'ascend');
    pos = pos(1:n);
end
    
end