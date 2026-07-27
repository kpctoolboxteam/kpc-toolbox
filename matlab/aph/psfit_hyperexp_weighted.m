function [MAP,Eapx,x] = psfit_hyperexp_weighted(E,n)
w = (E(end).^(1/length(E))).^-(1:length(E));
%% objective function weights
Eunscaled = E;
Escale = [];
for k = 1:length(E)
    Escale(k) = E(1).^k;
end
E = (E(:) ./ Escale(:))';
w = (log(E(end)).^(1/length(E))).^-(1:length(E));
%% initial point
x0 = rand(1,2*n); x0(1:n)=x0(1:n)/sum(x0(1:n));
%% optimization
global MAXITER;
global MAXCHECKITER;
global T0;
MAXITER=100;
MAXCHECKITER=1000;
TOL=1e-8;
EPSTOL=100*TOL;
options=optimset( 'Display','off', 'LargeScale','off','MaxIter',MAXITER, 'MaxFunEvals',1e10, ...
    'MaxSQPIter',500, 'TolCon',TOL, 'Algorithm', 'interior-point', 'OutputFcn', @psfit_hyperexp_weighted_outfun);
T0 = tic; % needed for outfun

%% optimization program
[x, f]=fmincon(@objfun,x0,[],[],[],[],x0*0+EPSTOL,[],@nnlcon,options);
[alpha,T] = topar(x);
MAP = {T,-T*ones(n,1)*alpha};
MAP = map_scale(MAP, Eunscaled(1));
MAP = map_normalize(MAP);
return

    function [alpha,T] = topar(x)
        alpha = x(1:n);
        T = diag(-1./(x((n+1):2*n)));
    end

    function [c,ceq] = nnlcon(x)
        [alpha,T] = topar(x);
        for j=1:3
            Eapx(j)=factorial(j)*alpha*(-inv(T)^j)*ones(n,1);
        end
        c = - alpha;
        ceq(1:3) = (w(1:3) .* abs( log(E(1:3)) - log(Eapx(1:3))))';
        ceq(end+1:end+n) = sum(alpha) - 1;
    end

    function f = objfun(x)
        [alpha,T] = topar(x);
        for j=1:length(E)
            Eapx(j)=factorial(j)*alpha*(-inv(T)^j)*ones(n,1);
        end
        f = w * abs( log(E) - log(Eapx) )';
    end

end