function [pi,t]=ctmc_transient(Q,pi0,t0,t1)
if nargin==2
    t1=pi0;
    t0=0;
    pi0=ones(1,length(Q));pi0=pi0/sum(pi0);
end
if nargin==3
    t1=t0;
    t0=0;
end
[t,pi]=ode23(@ctmc_transientode,[t0,t1],pi0);
%%[t,pi]=ode45(@ctmc_transientode,[t0,t1],pi0); % standard order 4-5
%[t,pi]=ode113(@ctmc_transientode,[t0,t1],pi0);

    function dpidt=ctmc_transientode(t,pi)
        pi=pi(:)';
        dpidt=pi*Q;
        dpidt=dpidt(:);
    end

end