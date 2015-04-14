function [PROB,Q]=ctmc_solve(Q,options)
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

if length(Q) > 6000
    warning('the order of Q is greater than 6000, i.e., %d elements. Press key to continue.',length(Q));
    pause;
end    

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
if nargin>1 & strcmpi(options,'plot')
    labels=cell(1,length(Qsav));
    descr=cell(1,length(Qsav));
    for i=1:length(Qsav)
        labels{i}=sprintf('%d,%0.3f',i,PROB(i));
        descr{i}=sprintf('');
        for j=1:length(Qsav)
            if Qsav(i,j)>0
            descr{i}=sprintf('%s%d,%0.2f\n',descr{i},j,Qsav(i,j));
            end
        end
    end
    graphlayout('adjMatrix',Qsav-diag(diag(Qsav)),'nodeLabels',labels,'nodeDescriptions',descr,'edgeLabels',labels,'currentLayout',treelayout);
end
