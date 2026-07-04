function [PH,scoreGamma,eflag,x,PHj] = kpcqbd_exact_fit(J,K,simprobs,rho)

gammaAllParams = 2;

warning off
opts = optimoptions(@fmincon,'Algorithm','interior-point', ...
    'MaxIter',1e3, ...
    'MaxFunEvals',1e10,  ...
    'MaxSQPIter',50,  ...
    'TolCon',1e-6,  ...
    'TolFun',1e-6,  ...
    'Display','off', ...
    'OutputFcn',@(x,optimValues,state) kpcfit_ph_outfun(x,optimValues,state));


fstPhUb = [Inf;Inf];

fstPhLb = [0;0];

A=[];
for j=1:J-1
    A = [A;zeros(1,gammaAllParams+(j-1)*(K-1)),ones(1,(K-1)),zeros(1,(J-1)*(2*K-1)-(j*(K-1)))];
end

x0 = genInitPoint(K,J);

problem = createOptimProblem('fmincon','x0',x0,...
        'objective',@objfun,'lb',[fstPhLb;zeros(length(x0)-2,1)],'ub',[fstPhUb;ones((K-1)*(J-1),1);Inf*ones(K*(J-1),1)], ...   
        'Aineq',A,'bineq',ones(J-1,1),...
        'options',opts);
[x,scoreGamma,eflag] = fmincon(problem);

[i,theta] = getGamma(x);
probs = reshape(x(gammaAllParams+1:gammaAllParams+(J-1)*(K-1)),K-1,J-1)';
lambdas = reshape(x(gammaAllParams+(J-1)*(K-1)+1:end),K,J-1)';
categparams = [probs,lambdas];

PHj{1} = fitAPHToGamma(i,theta);

for j=1:J-1
    p = categparams(j,1:K-1);
    pfull = [p,1-sum(p)];
    lambdas = 1./categparams(j,K:end);
    PHj{j+1} = {diag(-lambdas),lambdas'*pfull};
end

KPC_PH = PHj{1};

for i=2:J
    KPC_PH = map_kpc(KPC_PH,PHj{i});
end

PH{1} = KPC_PH;
PH{2} = scoreGamma;
PH{3} = x;


    function f=objfun(x)
       
        [k,theta] = getGamma(x);
        probs = reshape(x(gammaAllParams+1:gammaAllParams+(J-1)*(K-1)),K-1,J-1)';
        rates = reshape(x(gammaAllParams+(J-1)*(K-1)+1:end),K,J-1)';
        categparams = [probs,rates];

        currphlist{1} = fitAPHToGamma(k,theta);

        if ~map_isfeasible(currphlist{1})
            f = 1e10;
            return;
        end

        kpclambda = map_lambda(currphlist{1});

        for j=1:J-1
            p = categparams(j,1:K-1);
            pfull = [p,1-sum(p)];
            lambda_curr = 1./categparams(j,K:end);
            currphlist{j+1} = {diag(-lambda_curr),lambda_curr'*pfull};

            if sum(isnan(lambda_curr)) + sum(isinf(lambda_curr)) > 0 || sum(p) > 1
                f = 1e10;
                return;
            end
            kpclambda = kpclambda*map_lambda(currphlist{j+1});
        end

%         f=steadyStateError_exact(currphlist,simprobs,rho);
        f=steadyStateError_exact(currphlist,simprobs,rho*kpclambda);
        
    end

    function x0 = genInitPoint(phsize,noOfphs)
        pi0 = rand(1,(phsize-1)*(noOfphs-1));
        lambdaparams = rand(1,phsize*(noOfphs-1))*10;
        fstPhPars = [rand*10,rand];

        x0 = [fstPhPars,pi0,lambdaparams]';
    end


    function [k,theta] = getGamma(x)
        k = x(1);       
        theta = x(2);
    end

    function ph = fitAPHToGamma(k,theta)
        ph=aph_fit(k*theta,theta^2*gamma(k+2)/gamma(k),theta^3*gamma(k+3)/gamma(k),15);
    end
        
end