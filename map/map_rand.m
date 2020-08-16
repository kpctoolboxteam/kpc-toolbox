function [MAP]=map_rand(K)
if nargin<1
    K=2;
end
%D1=abs(sprandn(K,K,0.5));
D1=rand(K,K);
D0=rand(K,K);

MAP=cell(1,2);
MAP{1}=D0;
MAP{2}=D1;
MAP=map_normalize(MAP);
end


