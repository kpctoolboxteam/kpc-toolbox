function [C,bC]=trace_iat2bins(S,scale)
% computes the counts C of S in each bin with timescale "scale"
% binC(i) gives the bin of membership element i
n = length(S);
CS = cumsum(S);
bins = ceil((CS(end)-CS(1))/scale);
C = zeros(1,bins);
bC=[];
cur = 1;
last = 0;
for i=1:(bins+1)
    if cur == n
        break;
    end
    while CS(cur+1) <= i*scale
        cur = cur + 1;
        if cur == n
            break;
        end
    end
    C(i) = cur - last;
    bC(end+1:end+(cur - last)) = i;
    last = cur;
end
end