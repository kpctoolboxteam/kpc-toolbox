function M=det_acf(DET,kset)
P=det_embedded(DET);
n=length(DET{1});
M=[];
for k=kset
    for i=1:n
        for j=1:n
            K(i,j)=(-DET{1}(i,i))^-1*P(i,j);
        end
    end
    M(end+1)=(dtmc_solve(P)*K*P^(k-1)*K*e(n)-det_moment(DET,1)^2)/(det_moment(DET,2)-det_moment(DET,1)^2);
end
end