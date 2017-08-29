function Q=makeinfgen(Q)
for i=1:length(Q)
    Q(i,i)=-(sum(Q(i,:))-Q(i,i));
end
end