function Pnt=map_pntiter(MAP,na,t,M)
% Pnt=map_pntiter(MAP,n,t) - probability of having n arrivals within an interval of length t
% Neuts and Li, MAM1
if ~exist('M','var')
    M=ceil(log2(t*100/map_mean(MAP)));
end
if M<0
    [Pnt,P]=map_pntbisect(MAP,na,t);
else
    [Pnt,P]=map_pntbisect(MAP,na,t/2^M);
    for i=1:M
        Pold=P;
        for n=0:na
            P{n+1}=0*P{n+1};
            for j=0:n
                P{n+1}=P{n+1}+Pold{j+1}*Pold{n-j+1};
            end
        end
    end
    Pnt=P{na+1};
end
end

function [Pnt,P]=map_pntbisect(MAP,na,t)
% Pnt=map_pntiter(MAP,n,t) - probability of having n arrivals within an interval of length t
% Neuts and Li, MAM1
D0=MAP{1};
D1=MAP{2};
tau=max(diag(-D0));
N=findN(tau,t);
V=cell(na+1,N+1);
P=cell(na+1,1);
I=eye(size(D0));
K=D0/tau+I;
K1=D1/tau;
% initialize
for n=0:na
    P{n+1}=0*I;
    for k=0:N
        V{n+1,k+1}=0*I;
    end
end
% algorithm
V{0+1,0+1}=I;
P{0+1}=V{0+1,0+1}*br(tau,t,0);
for n=1:na
    V{n+1,0+1}=0*I;
    P{n+1}=0*I;
    for k=1:N
        V{n+1,k+1}=V{n+1,k-1+1}*K+V{n-1+1,k-1+1}*K1;
        P{n+1}=P{n+1}+V{n+1,k+1}*br(tau,t,n);
    end
end
Pnt=P{na+1};
end
function N=findN(tau,t)
epsilon=eps;
Nmax=100;
for N=1:Nmax
    S=0;
    for n=(N+1):Nmax
        S=S+br(tau,t,n);
    end
    if S<epsilon
        break;
    end
end
end
function b=br(tau,t,r)
b=exp(-tau*t)*(tau*t)^r/factorial(r);
end


