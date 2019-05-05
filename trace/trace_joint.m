function JM=trace_joint(S,lag,order)
% trace_joint - compute for a trace the joint moments E[X^{k_1}_{i} X^{k_2}_{i+j} ...]
% JM=trace_joint(S,lag,order)

lag=sort(cumsum(lag));
K=length(lag);
lag=lag-lag(1)*ones(1,K);

JMv = S(1:length(S)-max(lag)).^order(1);
for i = 2:length(order)
    JMv = JMv.*(S(lag(i)+1:length(S)-max(lag)+lag(i)).^order(i));
end
JM=mean(JMv);
end
