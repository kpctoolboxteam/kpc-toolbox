function [MAP,Eapx] = psfit_hyperexp(E,n)
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
    'MaxSQPIter',500, 'TolCon',TOL, 'Algorithm', 'interior-point', 'OutputFcn', psfit_hyperexp_outfun);
T0 = tic; % needed for outfun

%% optimization program
[x, f]=fmincon(@objfun,x0,[],[],[],[],x0*0+EPSTOL,[],@nnlcon,options);
[alpha,T] = topar(x);
MAP = {T,-T*ones(n,1)*alpha};
    function [alpha,T] = topar(x)
        alpha = x(1:n);
        T = diag(-1./(x((n+1):2*n)));
    end

    function [c,ceq] = nnlcon(x)
        [alpha,T] = topar(x);
        c = - alpha;
        ceq = sum(alpha) - 1;
    end

    function f = objfun(x)
        [alpha,T] = topar(x);
        for j=1:length(E)
            Eapx(j)=factorial(j)*alpha*(-inv(T)^j)*ones(n,1);
        end
        f = norm(log(E)-log(Eapx),1);
    end


end