function n=e(t)
if nargin<1
    n=[1;1];
else
    if iscell(t)
        n=ones(length(t{1}),1);
    else
    n=ones(t,1);    
    end
end