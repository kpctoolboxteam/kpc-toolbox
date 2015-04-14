function [alpha,T,PHR]=map2ph(MAPIN)
% [MAPOUT,A,pi]=map2ph(MAPIN) - Returns a PH distribution
%
%  Input:
%  MAPIN: a MAP in the form of {D0,D1}
%
%  Output:
%  alpha: entry probability vector of the PH distribution
%  T: subgenerator of the PH distribution
%  PHR: PH-renewal process with distribution (pi,A) in (D0,D1) notation

PHR=MAPIN;
T=MAPIN{1};
alpha=dtmc_solve(-inv(MAPIN{1})*MAPIN{2});
PHR{2}=MAPIN{2}*ones((length(MAPIN{2})),1)*pi;
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
