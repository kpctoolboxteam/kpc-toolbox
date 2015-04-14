function [lambda,prob]=ph2hyper(PH)
if (norm(PH{1}-diag(diag(PH{1})))>1e-10)
    error('The PH distribution is not hyper-exponential')
end
lambda=diag(-PH{1})';
prob=dtmc_solve(inv(-PH{1})*PH{2});
end

function PROB=ctmc_solve(Q,options)
% PROB=ctmc_solve(P) - Equilibrium distribution of the continuous-time
% Markov chain
%
%  Input:
%  P: infinitesimal generator of the continuous-time Markov chain
%
%  Output:
%  PROB: equilibrium distribution of the continuous-time Markov chain
%
%  Examples:
%  - ctmc_solve([-0.5,0.5;0.2,-0.2])
%
% MAP Queueing Networks Toolbox
% Version 1.0 	 15-Apr-2008
Qsav=Q;
zerocol=find(sum(abs(Q))==0);
if length(zerocol)>1
    warning('ctmc_solve: the infinitesimal generator is reducible');
    b=zeros(size(Q,1),1);
elseif length(zerocol)==0
    zerocol=size(Q,1);
end
b=zeros(size(Q,1),1); b(zerocol)=1;
Q(:,zerocol)=1; % add normalization condition
p=(Q')\b; % compute steady state probability vector
%p=bicgstab(Q',b); % compute steady state probability vector
PROB=p'; % p is a row vector
end

function PROB=dtmc_solve(P)
% PROB=dtmc_solve(P) - Equilibrium distribution of a discrete-time Markov
% chain
%
%  Input:
%  P: stochastic transition matrix of the discrete-time Markov chain
%
%  Output:
%  PROB: equilibrium distribution of the discrete-time Markov chain
%
%  Examples:
%  - dtmc_solve([0.5,0.5;0.2,0.8])
%
% MAP Queueing Networks Toolbox
% Version 1.0 	 15-Apr-2008
PROB=ctmc_solve(P-eye(size(P)));
end