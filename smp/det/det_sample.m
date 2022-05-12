function [SAMPLES,LAST,FIRST]=det_sample(DET,nSamples,initState,seed)

if nargin<3
%    pi=map_piq(DET); % time stationary initialization
    pi=map_pie(DET); % interval stationary initialization
    x=cumsum(pi);
    r=rand();
    initState=min(find(r<=x));
elseif length(initState)>1
    pi=initState;
    x=cumsum(pi);
    r=rand();
    initState=min(find(r<=x));    
end
RUNS=floor(nSamples/10000);
LS=[initState];
SAMPLES=[];
LAST=[];
FIRST=[];
for i=1:RUNS
    [S,L,F]=sub_map_sample(DET,10000,LS);
    SAMPLES(end+1:end+length(S),1)=S(:);
    LAST(end+1:end+length(L),1)=L(:);
    FIRST(end+1:end+length(L),1)=F(:);
end
[S,L,F]=sub_map_sample(DET,mod(nSamples,10000),LS);
SAMPLES(end+1:end+length(S),1)=S(:);
LAST(end+1:end+length(L),1)=L(:);
FIRST(end+1:end+length(L),1)=F(:);
end

function [SAMPLES,LAST,FIRST]=sub_map_sample(DET,nSamples,initState)
LAST=[];
nStates=length(DET{1});
for b=0:1
    for i=1:nStates
        for j=1:nStates
            p(i,b*nStates+j)=DET{b+1}(i,j)/abs(DET{1}(i,i));
        end
    end
end
for i=1:nStates
    p(i,i)=0;
end
cdf=0*p;
for i=1:nStates
    cdf(i,:)=cumsum(p(i,:));
end
cdf=abs(cdf);
curState=initState;
visits=zeros(nSamples,20);
maxpathlen=size(visits,2);
for i=1:nSamples
    arrival=0;
    last=2;
    visits(i,1)=curState;
    while ~arrival
        destState=find(cdf(curState,:)>=rand,1,'first');
        if destState>nStates
            arrival=1;
            destState=destState-nStates;
            curState=destState;
        else
            visits(i,last)=destState;
            curState=destState;
            last=last+1;
        end
        if last>maxpathlen
            maxpathlen=maxpathlen+5;
            visits(:,end+5)=0;
        end
    end
end
holdTimes=-1./diag(DET{1});
H=visits;
for i=1:nSamples
    LAST(i,1)=visits(i,max(find(visits(i,:))));
end
for i=1:nStates
    H(find(visits==i))=holdTimes(i);
end

SAMPLES=0*H(:,1);
for i=1:size(H,2)
    SAMPLES=SAMPLES+(H(:,i));
end
FIRST=visits(:,1);

end

