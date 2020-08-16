function MAP=map_sum(MAP,n)
for i=1:n
    MAPs{i}=MAP;
end
ns=[];
for i=1:n
    ns(i)=length(MAPs{i}{1});
end
D0=zeros(sum(ns));
D1=D0;

curpos=0;
for i=1:n
    D0((curpos+1):(curpos+ns(i)),(curpos+1):(curpos+ns(i))) = MAPs{i}{1};
    if i<n
        D0((curpos+1):(curpos+ns(i)),(curpos+ns(i)+1):(curpos+ns(i)+ns(i+1))) = MAPs{i}{2};
    else
        D1((curpos+1):(curpos+ns(i)),1:ns(1)) = MAPs{i}{2};
    end
    curpos = curpos + ns(i);
end
MAP={D0,D1};
end