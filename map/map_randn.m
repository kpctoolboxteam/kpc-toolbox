function [MAP]=map_randn(K,MU,SIGMA)
D0=abs(normrnd(MU(1),SIGMA(1),K,K));
D1=abs(normrnd(MU(1),SIGMA(1),K,K));
MAP=cell(1,2);
MAP{1}=D0;
MAP{2}=D1;
MAP=map_normalize(MAP);
end