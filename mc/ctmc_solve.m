function [p,Q]=ctmc_solve(Q,options)
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

if length(Q) > 6000 && ~options.force
    warning('the order of Q is greater than 6000, i.e., %d elements. Press key to continue.',length(Q));
    pause;
end

if size(Q)==1
    p = 1;
    return
end
n = length(Q);
zerocol=find(sum(abs(Q),1)==0);
nnzcol = setdiff(1:n, zerocol);
p = zeros(1,n);
b = zeros(n,1);
if length(zerocol)>=1
    warning('ctmc_solve: the infinitesimal generator is reducible');
end
Qnnz = Q(nnzcol, nnzcol);
bnnz = b(nnzcol);
Qnnz(:,end) = 1;
bnnz(end) = 1;
p(nnzcol)=Qnnz'\ bnnz;
end