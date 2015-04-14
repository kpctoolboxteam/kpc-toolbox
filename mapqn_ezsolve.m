function [XN,QN,UN,PROB,MAPQN]=mapqn_ezsolve(MAPS,N,P,ALPHA)
% [XN,QN,UN,PROB,MAPQN]=mapqn_ezsolve(MAP,N,P,alpha) - Define and solve
% using the underlying Markov chain a closed MAP queueing network with
% possibly load-dependent rates and arbitrary state-independent routing.
% Service is assumed always to be first-come first-served.
%
%  Input:
%  MAPS: a set of M MAPs describing the service process of each of the M
%  queues in the form {{D0_1,D1_1};{D0_2,D1_2};...;{D0_M,D1_M}}
%  N: model population, i.e., number of cyclic jobs
%  P: routing matrix, i.e., P(i,j) is the probability that a job departing
%  from i joins queue j (DEFAULT: cyclic network, all queues in series)
%  ALPHA: load-dependent service rate, where ALPHA(i,n) is the scaling
%  factor of queue i service rate when there are n jobs at i. Note that a
%  delay is defined by setting its ALPHA(i,n)=n. (DEFAULT: all ALPHAs equal
%  to one, load-independent model)
%
%  Output:
%  XN : XN is the mean throughput compute at the output of queue 1
%  QN : QN(i,h) is the queue-length of queue i when its service
%  process is in phase h
%  UN : UN(i,h) is the utilization of queue i when its service
%  process is in phase h
%  PROB : equilibrium probabilities of the network (see below)
%  MAPQN : data structure summarizing the network (see below)
%
%  MAPQN Data Structure:
%  MAPQN.Q : infinitesimal generator of the MAP queueing network
%  MAPQN.SS : all feasible distribution of jobs across the network
%  MAPQN.KK : all feasible combination of active phases across the service
%  processes
%  MAPQN.K : number of phases for each service process
%
%  PROB: equilibrium probabilities of states (NVEC, KVEC), where NVEC(i) is
%  the number of jobs in queue i and KVEC(i) is the current phase of its
%  service process. NVEC must be a row of MAPQN.SS, while KVEC must be a
%  row of MAPQN.KK. The probability p of state (NVEC, KVEC) is retrieved as
%
%  p=PROB(mapqn_hashstate(MAPQN,NVEC,KVEC))
%
%  Marginal probability can be computed immediately from these values.
%
%  Examples:
%  - mapqn_ezsolve({map_exponential(1);map_exponential(1)},10), throughput
%  of a product-form model with two balanced exponential servers in series 
%  and ten jobs
%  - mapqn_ezsolve({map_exponential(1);map_exponential(1)},10,[.1,.9;1,0]), 
%  same as above, but jobs self-loop at queue 1 with probability 0.1
%  - mapqn_ezsolve({map_exponential(1);map_exponential(1)},10,[.1,.9;1,0],
%  [ones(1,10);1:10]), same as above, but queue two is now a delay server
%  - mapqn_ezsolve({map_exponential(1);map_hyperexp(1,2)},10,[.1,.9;1,0],
%  [ones(1,10);1:10]), same as above, but queue two has now
%  hyperexponential service
%  - [XN,QN,UN,p,MAPQN]=mapqn_ezsolve({map_erlang(1,2);map_erlang(1,2)},4), 
%   -sum(QN(1,:)) is the mean queue-length at the exponential server,
%   -sum(UN(1,:)) is the utilization at the exponential server,
%   -p(mapqn_hashstate(MAPQN,[3,1],[2,2])) is the probability of the state
%    where queue 1 has 3 enqueued jobs and both erlang processes are in 
%    phase 2
%
%  NOTE: in the load-dependent case the first station must be constant rate
%
% MAP Queueing Networks Toolbox
% Version 1.0 	 15-Apr-2008
M=size(MAPS,1);
if nargin<4
    ALPHA=ones(M,N);
end
if nargin<3
    P = circul(e(M));
elseif isempty(P)
    P = circul(e(M));
end
MAPQN=mapqn_infgenrev(MAPS,N,P,ALPHA);
[XN,QN,UN,PROB]=mapqn_solverev(MAPQN,MAPS,ALPHA);
end
function MAPQN=mapqn_infgenrev(MAP,N,P,alpha)
% MAPQN_INFGENREV - Infinitesimal generator of closed MAP queueing network
% [Q,SS,KK,P]=mapqn_infgen(MAP,N,P,alpha)
% [XN,QN,UN,p,pmarg,QNb,pmargb]=mapqn_solve(Q,SS,muM,MAP,N)
% -- Input
% MAP     : cell array such that MAP{i} -> {D0,D1} is the MAP service process of queue i
% N       : network population (multiprogramming level)
% P       : routing matrix (alternative parameters: 'balanced' and 'cyclic')
% alpha   : load-dependent rate factor: qijhk -> qijhk*alpha(i,nvec(i)) (e.g., delay: alpha=1:N)
% --- Output
% Q       : infinitesimal generator matrix
% SS      : population state space

M=size(MAP,1);
%% compute routing matrix

%%
SS=multichoose(M,N); % population state space
K(1)=length(MAP{1}{1});
KK=1:K(1); KK=KK'; % phase state space
for i=2:M
    K(i)=length(MAP{i}{1}); % MAP order
    KKi=[];
    for k=1:K(i)
        KKi=[KKi;KK,k*ones(size(KK,1),1)];
    end
    KK=KKi;
end
clear KKi;

sz=size(SS,1)*size(KK,1);
Q=sparse(sz,sz);
for nveci=1:size(SS,1)
    nvec=SS(nveci,:);
    for kveci=1:size(KK,1)
        kvec=KK(kveci,:);
        row=(nveci-1)*size(KK,1)+kveci;
        % for all states (nvec,kvec) add q_{i,j}^{h,k} transitions
        for i=1:M % for all queue i, busy in phase h
            if nvec(i)>=1
                h=kvec(i); % h is the current phase
                for j=1:M % for all destinations
                    sspos=matchrow(SS,nvec+e(j)'-e(i)');
                    for k=1:K(i)
                        if i==j
                            qijhk=MAP{i}{1}(h,k)+P(i,i)*MAP{i}{2}(h,k);
                        else
                            qijhk=P(i,j)*MAP{i}{2}(h,k);
                        end
                        kkpos=matchrow(KK,kvec+(k-h)*e(i)');
                        col=(sspos-1)*size(KK,1)+kkpos;
                        Q(row,col)=Q(row,col)+qijhk*alpha(i,nvec(i));
                    end
                end
            end
        end
    end
end
Q=makeinfgen(Q);
MAPQN=struct('Q',Q,'SS',SS,'KK',KK,'K',K);
    function v=e(t)
        v=zeros(M,1); v(t)=1;
    end
end
function [XN,QN,UN,PROB]=mapqn_solverev(MAPQN,MAP,alpha)
M=length(MAP);
K=MAPQN.K;
Q=MAPQN.Q;
SS=MAPQN.SS;
KK=MAPQN.KK;
%
PROB=ctmc_solve(Q); % solve model
%PROB=ctmc_courtois(Q,{1:2:length(Q);2:2:length(Q)})
UN=zeros(M,max(K));
QN=zeros(M,max(K));
XN=0;
S1=0;
for nveci=1:size(SS,1)
    nvec=SS(nveci,:);
    for kveci=1:size(KK,1)
        kvec=KK(kveci,:);
        row=hash(nvec,kvec);
        for i=1:M
            h=kvec(i);
            if nvec(i)>0
                UN(i,h)=UN(i,h)+PROB(row);
                QN(i,h)=QN(i,h)+nvec(i)*PROB(row);
                if i==1
                    S1=S1+(map_mean(MAP{1})/alpha(1,nvec(1)))*PROB(row);
                end
            end
        end
    end
end
S1=S1/sum(UN(1,:)); % normalize PROB(row) to a conditional probability
XN=sum(UN(1,:))/S1;

    function col=hash(nvec,kvec)
        sspos=matchrow(SS,nvec);
        kkpos=matchrow(KK,kvec);
        col=(sspos-1)*size(KK,1)+kkpos;
    end
end
function Q=makeinfgen(Q)
for i=1:length(Q)
    Q(i,i)=-(sum(Q(i,:))-Q(i,i));
end
end
function [v] = multichoose(n,k,opt)
v=[];
if n==1
    v=k;
    return
elseif k==0
    v=zeros(1,n);
else
    last=0;
    for i=0:k
        w=multichoose(n-1,k-i);
        for j=1:size(w,1)
            v(end+1,:)=[i w(j,:)];
        end %for
    end %for
end %if

if nargin==3 & strcmpi(opt,'plot')
    G=zeros(size(v,1));
    for i=1:size(v,1)
        newrow=v(i,:);
        for j=1:size(v,2)
            if v(i,j)>0
                newrow(j)=newrow(j)-1;
                for k=1:size(v,2)
                    newrow(k)=newrow(k)+1;
                    d=matchrow(v,newrow);
                    if k==j+1 | (k==1 & j==size(v,2))
                        G(i,d)=1;
                    end
                    newrow(k)=newrow(k)-1;
                end
                newrow(j)=newrow(j)+1;
            end
        end
    end
    G
    labels=cell(1,size(v,1));
    for i=1:size(v,1)
        str='';
        for j=1:size(v,2)
            str=sprintf('%s%d',str,v(i,j));
        end
        labels{i}=str;
    end
    labels{1}
    labels{2}
end
end
function feasiblerows = matchrow(Matrix, row)
% feasiblerows = matchrow(Matrix, row)
feasiblerows=1:size(Matrix,1);
for col=1:length(row)
    A=Matrix(feasiblerows,col)==row(col);
    t=find(A);
    feasiblerows=feasiblerows(t);
    if length(feasiblerows)==1 && sum(Matrix(feasiblerows,:)==row)==size(Matrix,2)
        return
    end
end
feasiblerows=-1;
end
function C=circul(c)
if length(c)==1
    if c==1
        C=1;
    else
        v=zeros(1,c);
        v(end)=1;
        C=circul(v);
    end
    return
end
n=length(c);
I=eye(n);
R=I(1:n,[2:n 1]);
C=zeros(n);
for t=0:n-1
    C=C+c(1+t)*R^t;
end

end
function n=e(t)
if nargin<1
    n=[1;1];
else
    if iscell(t)
        n=ones(length(t{1}),1);
    else
    n=ones(t,1);    
    end
end
end