function [SAMPLES,LAST,FIRST]=map_sample(MAP,nSamples,pi,seed)
% [SAMPLES,LAST]=map_sample(MAP,NUM,S0, SEED) - Generate a random sample
%  of inter-arrival times
%
%  Input:
%  MAP: a MAP in the form of {D0,D1}
%  NUM: number of samples to be generated
%  S0: phase where the MAP starts for generating the first sample (DEFAULT:
%  random). if it is a vector it is interpreted as the probability of
%  starting from a phase
%  SEED: seed used by the random number generator
%
%
%  Output:
%  SAMPLES: set of N inter-arrival time samples
%  LAST: vector of absorbing phases for each sample
%  FIRST: vector of entering phase for each sample
%
%  Examples:
%  - [SAMPLES,LAST]=map_sample(MAP,10) - sample of 10 inter-arrival times
%  - [SAMPLES,LAST,FIRST]=map_sample(MAP,10,2) - sample of 10 inter-arrival
%  times starting from phase 2
%
%  * WARNING: the script can  be memory consuming and quite slow for
%  samples sizes greater than N=10000 *
%

if size(MAP{1})==1
    SAMPLES=exprnd(map_mean(MAP),nSamples,1);
    LAST=SAMPLES*0+1;
    FIRST=SAMPLES*0+1;    
    return
end
if nargin<3
    %    pi=map_piq(MAP); % time stationary initialization
    pi=map_pie(MAP); % interval stationary initialization
end
if nargin<4
    seed=rand('seed');
end
if length(pi)==1
    pi=zeros(1,length(MAP{1}));
    pi(initState)=1;
end
rand('seed',seed);
r=rand();
initState=min(find(r<=cumsum(pi)));
RUNS=floor(nSamples/10000);
LS=[initState];
SAMPLES=[];
LAST=[];
FIRST=[];
for i=1:RUNS
    [S,L,F]=sub_map_sample(MAP,10000,LS);
    SAMPLES(end+1:end+length(S),1)=S(:);
    LAST(end+1:end+length(L),1)=L(:);
    FIRST(end+1:end+length(L),1)=F(:);
end
[S,L,F]=sub_map_sample(MAP,mod(nSamples,10000),LS);
SAMPLES(end+1:end+length(S),1)=S(:);
LAST(end+1:end+length(L),1)=L(:);
FIRST(end+1:end+length(L),1)=F(:);
end

function [SAMPLES,LAST,FIRST]=sub_map_sample(MAP,nSamples,initState)
LAST=[];
nStates=length(MAP{1});
for b=0:1
    for i=1:nStates
        for j=1:nStates
            p(i,b*nStates+j)=MAP{b+1}(i,j)/abs(MAP{1}(i,i));
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
holdTimes=-1./diag(MAP{1});
H=visits;
for i=1:nSamples
    LAST(i,1)=visits(i,max(find(visits(i,:))));
end
for i=1:nStates
    H(find(visits==i))=holdTimes(i);
end
SAMPLES=0*H(:,1);
for i=1:size(H,2)
    SAMPLES=SAMPLES+exprnd(H(:,i));
end
FIRST=visits(:,1);
end
