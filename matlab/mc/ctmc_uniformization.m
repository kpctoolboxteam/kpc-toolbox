function [pi,kmax]=ctmc_uniformization(pi0,Q,t,tol,maxiter)
if ~exist('tol','var')
    tol = 1e-12;
end
if ~exist('maxiter','var')
    maxiter = 100;
end
q=1.1*max(abs(diag(Q)));
Qs=speye(size(Q))+sparse(Q)/q;
k=0;
s=1;
r=1;
iter=0;
kmax=1;
while iter<maxiter
    iter=iter+1;
    k=k+1;
    r=r*(q*t)/k;
    s=s+r;
    if (1-exp(-q*t)*s)<=tol
        kmax=k;
        break;
    end
end

pi=pi0*(exp(-q*t));
P=pi0;
ri=exp(-q*t);
for j=1:kmax
    P=P*Qs;
    ri=ri*(q*t/j);
    pi=pi+ri*P;
end
end
