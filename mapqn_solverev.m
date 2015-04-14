function [XN,QN,UN,C,p,p1,p1c,p2,COV]=mapqn_solverev(MAPQN,MAP)
% MAPQN_SOLVEREV - equilibrium performance of MAP closed queueing networks
% [XN,QN,UN,C,p,p1,p2]=mapqn_solverev(MAPQN)
% -- Input
% -- Output
% XN     : XN is the throughput of station 1
% QN     : QN(i,h) is the queue-length of station i in phase h
% UN     : UN(i,h) is the utilization of station i in phase h
% C ? conditional queue length?
% p       : p are the equilibrium probabilities of Q
% p1      : p1 linear reduction - busy period
% p1c     : p1c linear reduction - idle period
% p2      : p2 quadratic reduction
% COV     : queue-length covariance, COV(j,j)=VAR(j)
M=MAPQN.M;
K=MAPQN.K;
Q=MAPQN.Q;
SS=MAPQN.SS;
KK=MAPQN.KK;
P=MAPQN.P;
N=MAPQN.N;
%
p=ctmc_solve(Q);
if issym(p(1))
    p1=sym(zeros(M,max(K),M,N+1,max(K)));
    p1c=sym(zeros(M,max(K),M,N+1,max(K)));
    if nargout>7
        p2=sym(zeros(M,N+1,max(K),M,N+1,max(K)));
    end
    UN=sym(zeros(M,max(K)));
    QN=sym(zeros(M,max(K)));
    COV=sym(zeros(M,M));
    C=sym(zeros(M,max(K),M));
else
    p1=zeros(M,max(K),M,N+1,max(K));
    p1c=zeros(M,max(K),M,N+1,max(K));
    if nargout>7
        p2=zeros(M,N+1,max(K),M,N+1,max(K));
    end
    UN=zeros(M,max(K));
    QN=zeros(M,max(K));
    COV=zeros(M,M);
    C=(zeros(M,max(K),M));
end
for nveci=1:size(SS,1)
    nvec=SS(nveci,:);
    for kveci=1:size(KK,1)
        kvec=KK(kveci,:);
        row=hash(nvec,kvec);
        for i=1:M
            h=kvec(i);
            if nvec(i)>0
                UN(i,h)=UN(i,h)+p(row);
                QN(i,h)=QN(i,h)+nvec(i)*p(row);
            end
            for j=1:M
                k=kvec(j);
                if nvec(j)>0 % if busy
                    p1(j,k,i,nvec(i)+1,h)=p1(j,k,i,nvec(i)+1,h)+p(row);
                else % if idle
                    p1c(j,k,i,nvec(i)+1,h)=p1c(j,k,i,nvec(i)+1,h)+p(row);
                end
                if nargout>7
                   p2(j,nvec(j)+1,k,i,nvec(i)+1,h)=p2(j,nvec(j)+1,k,i,nvec(i)+1,h)+p(row);
                end
                if nargout>8
                    COV(i,j)=COV(i,j)+nvec(i)*nvec(j)*p(row);
                end
            end
        end
    end
end
XN=sum(UN(1,:))/map_mean(MAP{1});
if nargout>3
    for j=1:M
        for k=1:K(j)
            for i=1:M
                C(j,k,i)=0;
                for h=1:K(i)
                    for ni=1:N
                        C(j,k,i)=C(j,k,i)+ni*p1(j,k,i,ni+1,h);
                    end
                end
            end
        end
    end
end
    function col=hash(nvec,kvec)
        sspos=matchrow(SS,nvec);
        kkpos=matchrow(KK,kvec);
        col=(sspos-1)*size(KK,1)+kkpos;
    end

end