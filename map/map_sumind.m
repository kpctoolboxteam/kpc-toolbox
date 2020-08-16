function MAP=map_sumind(MAPs)
n=length(MAPs); % MAPs to be summed
l=[];
for i=1:n
    l(i)=length(MAPs{i}{1});
end
D0=zeros(sum(l));
D1=D0;

curpos=0;
for i=1:n
   D0((curpos+1):(curpos+l(i)),(curpos+1):(curpos+l(i))) = MAPs{i}{1};       
   if i<n
       D0((curpos+1):(curpos+l(i)),(curpos+l(i)+1):(curpos+l(i)+l(i+1))) = MAPs{i}{2}*ones(l(i),1)*map_pie(MAPs{1+i});
   else       
       D1((curpos+1):(curpos+l(i)),1:l(1)) = MAPs{i}{2}*ones(l(i),1)*map_pie(MAPs{1});
   end
   curpos = curpos + l(i);
end
MAP={D0,D1};
end