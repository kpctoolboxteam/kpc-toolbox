function f=map_pdf(MAP,tset)
% MAP_PDF computes the probability density function (PDF) of a MAP at the point
% tset
% MAP:  Markovian arrival process to evaluate
% tset: time-point to evaluate the PDF of the MAP
%
% Copyright (c) 2012-2014, Imperial College London
% All rights reserved.


pi=map_pie(MAP);
f=[];
if issym(MAP{1})
    f=sym(f);
end
e=ones(length(MAP{2}),1);
if year(matlabRelease.Date)>2023 || strcmpi(matlabRelease.Release,"R2023b")
    f = expmv(MAP{1}',pi',tset)'*(-MAP{1})*e;
    f = f';
else
    for t=tset
        f(end+1)=pi*(expm(MAP{1}*t)*(-MAP{1}))*e;
    end
end
end

