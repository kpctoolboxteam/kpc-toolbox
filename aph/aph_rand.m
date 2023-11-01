function [MAP]=aph_rand(K)
if nargin<1
    K=2;
end
%D1=abs(sprandn(K,K,0.5));
D1=rand(K,K);
D0=rand(K,K);
for k=1:K
    D0(k,1:k-1)=0;
end
MAP=cell(1,2);
MAP{1}=D0;
MAP{2}=D1;
MAP=map_normalize(map_renewal(MAP));
end


