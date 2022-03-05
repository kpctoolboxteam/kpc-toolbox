function C=trace_iat2counts(S,scale)
% computes the counting process of S, i.e., the counts after "scale" units
% of time from an arrival
n = length(S);
CS = cumsum(S);
C = zeros(1,n-1);
for i=1:n-1
    cur = i;
    while CS(cur + 1) - CS(i) <= scale
        cur = cur + 1;
        if cur == n % when the window first hits the end of the trace we return
            C(i) = cur - i;            
            C=C(1:i);
            return;
        end
    end
    C(i) = cur - i;
end
end
