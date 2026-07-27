function [pi,t]=ctmc_transient(Q,pi0,t0,t1,useStiff,reltol)
if ~exist('useStiff','var')
    useStiff = false;
end
if ~exist('reltol','var')
    reltol = 1e-3;
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


odeoptions = odeset('RelTol', reltol);
if isinf(t1)
    nonZeroRates = abs(Q(Q~=0));
    nonZeroRates = nonZeroRates( nonZeroRates >0 );
    T = abs(100/min(nonZeroRates)); % solve ode until T = 100 events with the slowest rate
    if useStiff
        [t,pi]=ode15s(@ctmc_transientode, [t0,T], pi0, odeoptions);
    else
        [t,pi]=ode23(@ctmc_transientode, [t0,T], pi0, odeoptions);
    end
else
    if useStiff
        [t,pi]=ode15s(@ctmc_transientode, [t0,t1], pi0, odeoptions);
    else
        [t,pi]=ode23(@ctmc_transientode, [t0,t1], pi0, odeoptions);
    end
end
%[t,pi]=ode113(@ctmc_transientode,[t0,t1],pi0);

function dpidt=ctmc_transientode(t,pi)
pi=pi(:)';
dpidt=pi*Q;
dpidt=dpidt(:);
end

end