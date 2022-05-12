function M=det_moment(DET,kset)
n=length(DET{1});
M=[];
for k=kset
M(end+1)=dtmc_solve(det_embedded(DET))*inv(-DET{1})^k*e(n);
end
end