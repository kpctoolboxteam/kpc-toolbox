function [pi,t]=ctmc_transient(Q,pi0,t0,t1,useStiff)
if ~exist('useStiff','var')
useStiff = false;
end
if nargin==2
    t1=pi0;
    t0=0;
    pi0=ones(1,length(Q));pi0=pi0/sum(pi0);
end
if nargin==3
    t1=t0;
    t0=0;
end
if useStiff
    [t,pi]=ode15s(@ctmc_transientode,[t0,t1],pi0);
else
    [t,pi]=ode23(@ctmc_transientode,[t0,t1],pi0);
end
%[t,pi]=ode113(@ctmc_transientode,[t0,t1],pi0);

    function dpidt=ctmc_transientode(t,pi)
        pi=pi(:)';
        dpidt=pi*Q;
        dpidt=dpidt(:);
    end

end