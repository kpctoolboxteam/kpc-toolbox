function [MMPP]=mmpp_rand(K)
D0=rand(K,K);
D1=rand(K,K);
D1=diag(diag(D1));
MMPP=cell(1,2);
MMPP{1}=D0;
MMPP{2}=D1;
MMPP=map_normalize(MMPP);
end


