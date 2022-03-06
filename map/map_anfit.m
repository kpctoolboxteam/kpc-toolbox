function MAP=map_anfit(ls,rho,H,n,ds,SA,SAlags, iter_max, iter_tol)
% MAP_ANFIT - Andersen and Nielsen MAP Fitting Algorithm
% From A.T.Andersen and B.F.Nielsen, "A Markovian Approach for Modeling Packet Traffic with Long-Range Dependence", IEEE JSAC 16(5), 1998.
% ls - mean arrival rate
% rho - lag-1 autocorrelation of the counting process
% H - Hurst coefficient (counting process)
% n - number of time scales to be modeled
% ds - number of IPP to be used 
% SA - set of autocorrelation coefficients of interarrival times
% SAlags - lags used in SA
%
% EXAMPLE: 
% lstar=1; rho=0.8; H=0.92; n=5; ds=4; SA=[0.3,0.28]; SAlags=[1,2];
% MAP=map_anfit(lstar,rho,H,n,ds,SA,SAlags);

LSQFIT=0;
if nargin>5 % if SA and SAlags are given in input, then perform least square fit of SA
    LSQFIT=1;
    SA=SA(SAlags);
end
if nargin<6
        iter_max=100;
end       
if nargin<7
        iter_tol=1e-9;
end       

%%
beta=2-2*H;
d=ds;
k(2,1)=0.8;
%%
skip=0;
while 1
    if ~skip
        d0=0; phi(d)=1; i=1; a=10^(n/(d-1));
    end
    S=0;
    for j=0:(i-1)
        S=S+phi(d-j)^2*exp(1-a^(i-j));
    end
    D=a^(i*beta)-S;
    if D<0
        phi(d-i)=0;
        d0=d0+1;
        if ds>d-d0
            d=d+1;
            skip=0;
            continue;
        else
            i=i+1;
            if i==d
                break
            else
                skip=1;
                continue;
            end
        end
    else
        phi(d-i)=sqrt(D);
        i=i+1;
        if i==d
            break
        else
            skip=1;
            continue;
        end
    end
end
%%
for i=2:d
    k(2,i)=a^(1-i)*k(2,1);
end
if ~(k(2,1)<1 && rho<0.5)
    fprintf('warning: if necessary adjust k(2,1) and/or rho\n')
end
S=0;
for i=1:d
    kappa=k(2,i);
    e=exp(-kappa);
    S=S+phi(i)^2*kappa^(-2)*((1-e)^2-2*rho*(kappa-(1-e)));
end
eta=sqrt(4*rho*ls)/sqrt(S);
L=eta*sum(phi(1:d))/2;
if ls<L
    lP=0;
    for i=1:d
        c(1,i)=L^2/(ls^2+L^2)*k(2,i);
        c(2,i)=k(2,i)-c(1,i);
        l(i)=phi(i)*(ls^2+L^2)/(ls*sum(phi(1:d)));
    end
else
    lP=ls-L;
    for i=1:d
        c(2,i)=0.5*k(2,i);
        c(1,i)=c(2,i);
        l(i)=eta*phi(i);
    end
end


if ~LSQFIT
    %% No Least Square Fit, compose and return
    MAP=map_exponential(lP);
    for i=1:d
        D0=[0,c(1,i);c(2,i),0];
        D1=[l(i),0;0,0];
        IPP=map_normalize({D0,D1});
        MAP=map_super(MAP,IPP);
    end
    MAP=map_normalize(MAP);
    return
else
    %% Least Square Fit
    NSA=norm(SA,2);
    for i=1:d
        k(1,i)=(l(i)-0)^2*(c(1,i)*c(2,i))/(c(1,i)+c(2,i))^3;
        ls(i)=(c(2,i)*l(i))/(c(1,i)+c(2,i));
    end
    MAP0=[];
    r0=[];
    for i=1:d
        r0(end+1)=1+k(1,i)*k(2,i)/ls(i)^2;
    end
    [f,r]=minimize(@objfun,@nnlcon,r0, iter_max, iter_tol);
    MAP=MAP0; % return the last MAP found in the lsq algorithm
end
    function [c,ceq]=nnlcon(r)
        c=[];
        ceq=[];
        for i=1:d
            c(end+1)=-(ls(i)-sqrt(k(1,i)*k(2,i)/r(i)));
        end
    end
    function f=objfun(r)
        MAP0={[-lP],[lP]};
        for i=1:d
            D0=[0,k(2,i)*r(i)/(1+r(i));k(2,i)/(1+r(i)),0];
            D1=[ls(i)+sqrt(k(1,i)*k(2,i)*r(i)),0;0,ls(i)-sqrt(k(1,i)*k(2,i)/r(i))];
            IPP=map_normalize({D0,D1});
            MAP0=map_super(MAP0,IPP);
        end
        f=norm(map_acf(MAP0,SAlags)-SA',2);
    end

end


%% MINIMIZE encapsulates MATLAB's nonlinear optimization solver
% objfun  - objective function
% nnlcon  - nonlinear constraints
% x0      - initial solution
% MAXITER - maximum number of iterations
% TOL     - constraint feasibility tolerance
function [f,x]=minimize(objfun,nnlcon,x0,MAXITER,TOL)
warning off;
options = optimset('LargeScale','off','MaxIter',MAXITER, 'MaxFunEvals',1e10, 'MaxSQPIter',500, 'TolCon',TOL, 'Display','off');
EPSTOL=100*TOL;
xlb=EPSTOL+0*x0;
[x,f]=fmincon(objfun,x0,[],[],[],[],xlb,[],nnlcon,options);
end

%% MAP_NORMALIZE normalizes the diagonal element of D0 and prevents the solver to step into complex numbers
% MAP     - 2x1 cell, element 1 is the matrix D0, element 2 is D1
function [MAP]=map_normalize(MAP)
nStates=size(MAP{1},1);
nMAP=length(MAP);

for n=1:nStates
    MAP{1}(n,n)=0;
    for b=1:nMAP
        MAP{1}(n,n)=MAP{1}(n,n)-sum(MAP{b}(n,:));
    end
end
for i=1:nStates
    for j=1:nStates
        MAP{1}(i,j)=real(MAP{1}(i,j));
        MAP{2}(i,j)=real(MAP{2}(i,j));
    end
end
end

%% MAP_SUPER creates the superposed process MAPa \oplus MAPb
%% (oplus=Kronecker sum)
% MAPa     - 2x1 cell, element 1 is the matrix D0, element 2 is D1
% MAPb     - 2x1 cell, element 1 is the matrix D0, element 2 is D1
function MAP=map_super(MAPa,MAPb)
D0s=krons(MAPa{1},MAPb{1});
D1s=krons(MAPa{2},MAPb{2});
MAP={D0s,D1s};
MAP=map_normalize(MAP);
end

%% KRONS - Kronecker sum
function AB=krons(A,B)
AB=kron(A,eye(size(B)))+kron(eye(size(A)),B);
end

%% MAP_PROB - MAP equilibrium probabilities from D0 and D1
function p=map_prob(MAP)
Q=MAP{1}+MAP{2};
p=ctmc_solve(Q);
end

%% CTMC_SOLVE solve Markov process by global balance
function p=ctmc_solve(Q)
% Q - infinitesimal generator of the process
% p - steady state probabilities
zerocol=find(sum(abs(Q))==0);
if length(zerocol)>1
    warning('two cols of Q are composed by zeros only.');
    b=zeros(size(Q,1),1);
elseif length(zerocol)==0
    zerocol=size(Q,1);
end
b=zeros(size(Q,1),1); b(zerocol)=1;
Q(:,zerocol)=1; % add normalization condition
p=(Q')\b; % compute steady state probability vector
p=p'; % p is a row vector
end

%% MAP_ACF autocorrelation coefficents of inter-arrival times
% MAP     - 2x1 cell, element 1 is the matrix D0, element 2 is D1
% lagset  - set of lag coefficients
function acfun=map_acf(MAP,lagset)
p=map_prob(MAP); %steady state probabilities
lambda=p*MAP{2}*ones(size(MAP{2},2),1); % arrival rate
e=ones(size(MAP{1},1),1);
invD0=inv(-MAP{1});
acfun=[];
P=invD0*MAP{2};
for lag=lagset
    acfun(end+1)=lambda*p*(P^lag)*invD0*e;
end
acfun=(acfun-1)./map_scv(MAP);
end
