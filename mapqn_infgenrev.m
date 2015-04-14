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
if nargin<4
    if issym(MAP{1}{1})
        alpha=sym(ones(M,N));
    else
        alpha=ones(M,N);
    end
end
if nargin<3
    P = circul(e(M));
elseif isempty(P)
    P = circul(e(M));
end

%%
SS=multichoose(M,N); % population state space
%SS=sortbynnzpos(SS);
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
    if issym(MAP{1}{1})
    Q=sym(zeros(sz,sz));
    else
        Q=sparse(sz,sz);
    end
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
MAPQN=struct('Q',Q,'SS',SS,'KK',KK,'P',P,'K',K,'M',M,'N',N,'alpha',alpha);
    function v=e(t)
        v=zeros(M,1); v(t)=1;
    end
end


