function P=dtmc_makestochastic(P)
for i=1:size(P,1)
    if sum(P(i,:))>0
        P(i,:)=P(i,:)/sum(P(i,:));
        P(i,i)=1-(sum(P(i,:))-P(i,i));
    else
        P(i,:)=0;
        P(i,i)=1;
    end
end
end