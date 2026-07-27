function [resSCV,resG2,fobjAC]=kpcfit_sub_acfit(E,SA,SAlags,J,MaxIterAC, MaxRunsAC, MaxResAC, AnimateAC)
%% optimization parameters
MaxTimeAC = 10; % 10sec
MaxResAC=min([MaxResAC,MaxRunsAC]);


%% other variables
SCV=(E(2)-E(1)^2)/E(1)^2;
NSA=norm(SA,2); % normalization constant for objective function
SCVJ=SCV; % initialize nnlcon variable
%%
global MAXITER;
global MAXCHECKITER;
global T0;
global lgkx;
global stagnval;
global stagniter;
lgkx = [];
stagnval = 0;
stagniter = 0;

%% multi-start optimization
fset = [];   % record the objective function values for each param set
xparamset = {}; % record all results
xparam_best = [];
f_best = inf;
bestindex = 0;
for i = 1:MaxRunsAC
    x0=[(1+rand(J,1)); rand(J,1)]; x0=x0(:)';
    tstart = tic();
    optimoptions = optimset('Algorithm','interior-point', ...
        'LargeScale','off', ...
        'MaxIter',MaxIterAC, ...
        'MaxFunEvals',1e10,  ...
        'MaxSQPIter',50,  ...
        'TolCon',kpcfit_tol,  ...
        'Display','off', ...
        'OutputFcn',@(x,optimValues,state) kpcfit_sub_acfit_outfun(x,optimValues,state,MaxIterAC,MaxTimeAC,tstart,f_best));
    fprintf(1,'acfit: run %d of %d ',i,MaxRunsAC);
    tic;
    [x,f]=fmincon(@objfun,x0,[],[],[],[],0*x0+kpcfit_tol,[],@nnlcon,optimoptions);
    fprintf(1,'- objfun: %d (%3.3f sec) ',f,toc);
    xparamset{end+1} = x;
    fset(1,end+1) = f;
    if f < f_best
        fprintf(1,'**best**',toc);
        f_best = f;
        xparam_best = x;
    end
    fprintf(1,'\n',toc);
end
[vl,ind] = sort(fset,2,'ascend');

resG2 = {};
resSCV = {};
fobjAC = [];
for i = 1:MaxResAC
    [SCVj,G2j]=xtopar(xparamset{ind(i)});
    G2j(find(G2j>1-kpcfit_tol))=1-kpcfit_tol;
    resSCV{end+1} = SCVj;
    resG2{end+1} = G2j;
    fobjAC(i) = fset(1,ind(i));
end

%%
    function [SCVj,G2j]=xtopar(x)
        SCVj = x(1:J);
        G2j  = x((J+1):end);
    end

    function [c,ceq]=nnlcon(x)
        [SCVj,G2j]=xtopar(x);
        ceq=[];
        c=[];
        %% SCV CONSTRAINTS
        c(end+1) = (0.5-kpcfit_tol)-SCVj(1);    % SCVj(1) > 0.5
        for j=2:J
            c(end+1) = (1+kpcfit_tol) - SCVj(j);         % SCVj(j) > 1.1
        end
        %% G2 DECAY RATE CONSTRAINTS
        for j=1:J
            c(end+1) = G2j(j) - (1-kpcfit_tol);     % G2j(j) < 1
            c(end+1) = kpcfit_tol - G2j(j);       % G2j(j) > 0
        end
    end

    function f=objfun(x)
        [SCVj,G2j]=xtopar(x);
        [SCVJ,acfCoeff]=kpcfit_sub_eval_acfit(SCVj,G2j,SAlags);
        if AnimateAC
            loglog(SA); hold all; loglog(acfCoeff);  ylim([1e-6,1]);
            ylabel('autocorrelation \rho_k ');
            xlabel('lag k - [log]'); hold off;
            pause(0.001)
        end
        f=norm((SA-acfCoeff),1)/NSA + (SCVJ - SCV)^2/SCV^2;
    end

end