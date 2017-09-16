function [KPC_PH,score,x] = kpcfit_ph_search(E,J,options,x0,max_aph_order)
% Not for direct call, auxiliary function. Optimization-based search.
%% initialization
warning off
%optimoptions = optimset('Algorithm','active-set', ...
%optimoptions = optimset('Algorithm','trust-region-reflective', ...
optimoptions = optimset('Algorithm','interior-point', ...
    'LargeScale','off', ...
    'MaxIter',200, ...
    'MaxFunEvals',1e10,  ...
    'MaxSQPIter',50,  ...
    'TolCon',kpcfit_tol,  ...
    'Display','off', ...
    'OutputFcn',@(x,optimValues,state) kpcfit_ph_outfun(x,optimValues,state));

if nargin  < 5
    SCV = (E(2)-E(1)^2)/E(1)^2;
    max_aph_order = ceil(1/SCV);
end
global lgkx;
global stagnval;
global stagniter;
stagnval = 0;
stagniter = 0;
order = 2*ones(1,J); % Current version assumes all PH(2)s
%% x0 initial point
J = length(order); % number of PHs to compose
if nargin < 4 || length(x0) == 0
    x0 = zeros(sum(2*order-1),1); % each PH has 2n-1 degrees of freedom
end
w = (log(E(end)).^(1/length(E))).^-(1:length(E));
lgkx = x0; % last good known x

%reg = regress(log10(E(:)), [1:length(E)]'); % regression of growth trend of empirical moments
for j = 1:J
    PH{j} = map_feasblock(rand,1000*rand,-1,0); % hyperexponential
end
K = max([2*max(order)-1,length(E)]);
F = factorial(1:K);
logEtable = zeros(J,K);
for k = 1:K
    for j = 1:J
        logEtable(j,k) = log(map_moment(PH{j},k));
    end
end
if ~(nargin < 4 || length(x0) == 0)
    x0 = logEtable(:);
end
%%
[x,score] = fmincon(@objfun,x0,[],[],[],[],-200*ones(size(x0)),200*ones(size(x0)),@nnlcon,optimoptions); % 200 avoids exp(x) to overflow/underflow
logEtable = reshape(x,J,K);
[c,ceq]=nnlcon(x);
j=1;
if SCV<1
    PH{1} = map_scale(aph_fit(exp(logEtable(j,1)),exp(logEtable(j,2)),exp(logEtable(j,3)),max_aph_order),exp(logEtable(j,1)));
else
    PH{1}=map_feasblock(exp(logEtable(j,1)),exp(logEtable(j,2)),exp(logEtable(j,3)),0);
end
for j=2:J
    PH{j}=map_feasblock(exp(logEtable(j,1)),exp(logEtable(j,2)),exp(logEtable(j,3)),0);
end
KPC_PH = PH{1};
for j=2:J
    KPC_PH = map_kpc(KPC_PH,PH{j});
end
%keyboard
return

    function [c,ceq]=nnlcon(x)
        logEtable = reshape(x,J,K);
        logEcur = sum(logEtable,1) - (J-1)*log(F);
        c=[];
        
        j = 1;
        if SCV < 1
            j = 1;
            PH{1} = map_scale(aph_fit(exp(logEtable(j,1)),exp(logEtable(j,2)),exp(logEtable(j,3)),max_aph_order),exp(logEtable(j,1)));
            logEtable(j,1:K) = log(map_moment(PH{j},1:K));
            for j=2:J
                PH{j}=map_scale(map_feasblock(exp(logEtable(j,1)),exp(logEtable(j,2)),exp(logEtable(j,3)),0),exp(logEtable(j,1)));
                logEtable(j,1:K) = log(map_moment(PH{j},1:K));
            end
            
            % APH moment constraints
            % log e3 - log e1 -log e2 >= log e2 - 2*log e1 + log (n+1) -log n
            % log n2 >= log n - log (n-1)
            c(end+1) = -(logEtable(j,2) - log(max_aph_order) + log(max_aph_order-1));
            c(end+1) = -(logEtable(j,3) - 2*logEtable(j,2) +logEtable(j,1) - log(max_aph_order+1) + log(max_aph_order));
        else
            c(end+1) = -(logEtable(j,2) - 2*logEtable(j,1) - log(2)) - 10*kpcfit_tol; % SCV > 1
            c(end+1) = -(logEtable(j,3) - (log(3/2) + 2*logEtable(j,2) - logEtable(j,1))); % E3 > (3/2)*E2^2/E1
        end
        
        for j=2:J % for hyper-exponential PHs
            c(end+1) = -(logEtable(j,2) - 2*logEtable(j,1) - log(2)) - 10*kpcfit_tol; % SCV > 1
            c(end+1) = -(logEtable(j,3) - (log(3/2) + 2*logEtable(j,2) - logEtable(j,1))); % E3 > (3/2)*E2^2/E1
        end
        
        ceq = w(1:options.MinExactMom) .* ( log(E(1:options.MinExactMom)) - logEcur(1:options.MinExactMom)) ;
        
        if sum(isnan(c))>0
            c = 1e10*ones(size(c));
        end
        
        if sum(isnan(ceq))>0
            ceq = 1e10*ones(size(ceq));
        end
        
    end

    function f=objfun(x)
        logEtable = reshape(x,J,K);
        if isnan(x)
            f = 1e10;
            return
        end
        
        if SCV < 1
            j = 1;
            PH{1} = map_scale(aph_fit(exp(logEtable(j,1)),exp(logEtable(j,2)),exp(logEtable(j,3)),max_aph_order),exp(logEtable(j,1)));
            logEtable(j,1:K) = log(map_moment(PH{j},1:K));
            for j=2:J
                PH{j}=map_scale(map_feasblock(exp(logEtable(j,1)),exp(logEtable(j,2)),exp(logEtable(j,3)),0),exp(logEtable(j,1)));
                logEtable(j,1:K) = log(map_moment(PH{j},1:K));
            end
        else
            % determine moments using characteristic polynomial method
            for j=2:J
                m{j} = kpcfit_hyper_charpoly(exp(logEtable(j,:))./F,2)';
                if sum(isnan(m{j})) > 0 % if the result is NaN most likely we hit the exponential distribution
                    for k = 4:K % recursively compute moments of order 4 or more
                        idxs = (k-1):-1:(k-(2*order(j)-2));
                        Ej = map_moment({[-1],[1]},idxs);
                        logEtable(j,k) = log(-F(k)*m{j}(2:end)*(Ej./F(idxs))');
                    end
                else
                    for k = 4:K % recursively compute moments of order 4 or more
                        idxs = (k-1):-1:(k-(2*order(j)-2));
                        Ej = exp(logEtable(j,idxs));
                        if isnan(Ej) % reset this PH to exponential
                            Ej = map_moment({[-1],[1]},idxs);
                        end
                        logEtable(j,k) = log(-F(k)*m{j}(2:end)*(Ej./F(idxs))');
                    end
                end
            end
        end
        logEcur = sum(logEtable,1) - (J-1)*log(F);
        f= w * abs( log(E) - logEcur )';
        if isnan(f)
            f = 1e10;
        end
    end

end

function [MAP]=map_rand(K)
if nargin<1
    K=2;
end
%D1=abs(sprandn(K,K,0.5));
D1=rand(K,K);
D0=rand(K,K);

MAP=cell(1,2);
MAP{1}=D0;
MAP{2}=D1;
MAP=map_normalize(MAP);
end