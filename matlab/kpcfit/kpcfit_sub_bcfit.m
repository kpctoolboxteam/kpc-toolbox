function [E1j,E3j,f]=kpcfit_sub_bcfit(E,SCVj,G2j,BC,BCLags,MaxIterBC, MaxRunsBC)
%% optimization parameters
MaxTimeBC = 100;
NumMAPs=length(SCVj);
TOL=1e-9;
EPSTOL=10*TOL;
VERBOSE=0;


%% other variables
NBC=norm(BC,2); % normalizing constants for bicovariances

% truncate bicorrelations when it gets negative
tag=[];
for index=1:size(BCLags,1)
    if find(diff(BCLags(index,:))<0)
        tag(end+1,1)=index;
    end
end
BCLags(tag,:)=[];
BC(tag)=[];

E1j=E(1).^(1/NumMAPs)*ones(1,NumMAPs);
t=E(1)*rand;
fold=0;

xparamset = {}; % record all the possible x params
fset = [];   % record the objective function values for each param set
xparam_best = [];
f_best = inf;
for index = 1:NumMAPs
    E2j(index)=(1+SCVj(index))*E1j(index)^2;
    E3j(index)=(3/2+t)*E2j(index)^2/E1j(index);
end
x0=[E1j(:);E3j(:)]';
x0base=x0;


%fprintf(1,'bcfit: maximum number of solver iterations per run is MaxIterBC = %d\n', MaxIterBC);
for ind = 1:MaxRunsBC
    tstart = tic();
    if VERBOSE ~= 0
        fprintf(1,'bcfit: run %d of %d ',ind, MaxRunsBC);
    end
    tic;

    options = optimset('LargeScale','off', ...
        'MaxIter',MaxIterBC, ...
        'MaxFunEvals',1e10, ...
        'MaxSQPIter',500, ...
        'TolCon',TOL, ...
        'Display','off', ...
        'DiffMinChange',1e-10, ...
        'OutputFcn',@(x,optimValues,state) kpcfit_sub_bcfit_outfun(x,optimValues,state,MaxIterBC,MaxTimeBC,tstart,f_best));
    
    [x,f,exitflag,output]=fmincon(@objfun,x0,[],[],[],[],0*x0+EPSTOL,[],@nnlcon,options);
    
    if VERBOSE ~= 0
        fprintf(1,'(%3.3f sec, %d iter) f=%f',toc, output.iterations, f);
    end
    
    xparamset{end+1} = x;
    fset(1,end+1) = f;
    if f < f_best
        if VERBOSE ~= 0
            fprintf(1,'**best**',toc);
        end
        f_best = f;
        xparam_best = x;
    end
    
    if VERBOSE ~= 0
        fprintf(1,'\n',toc);
    end
   
   x0 = x0base .* [(0.25+1.75*rand(NumMAPs,1)); ones(NumMAPs,1)]';
   t = rand*E(1);
   for index = 1:NumMAPs
        E2j(index)=(1+SCVj(index))*x0(index)^2;
        x0(NumMAPs+index)=(3/2+t)*E2j(index)^2/x0(index);
    end
end

[E1j,E3j]=xtopar(xparam_best);
E1j(1)=E(1)/prod(E1j(2:end));
prod(E1j);
f = f_best;


    function [E1j,E3j]=xtopar(x)
        E1j=x(1:NumMAPs);
        E3j=x((NumMAPs+1):end);
        E1j(1)=E(1)/prod(E1j(2:end));
    end
    function [c,ceq]=nnlcon(x)
        [E1j,E3j]=xtopar(x);
        c=[];
        ceq=[];
        for j=2:NumMAPs
            c(end+1)=(2+EPSTOL)*E1j(j)^2-E2j(j);    % E2j(j) > 2*E1j(j)^2
            c(end+1)=(3/2+EPSTOL)*E2j(j)^2/E1j(j)-E3j(j);
            % E3j(j) >
            % 3*E2j(j)^2/(2*E1j(j))
        end
        if SCVj(1)>1
            % if SCV for the first MMPP2 is > 1
            % add the constraint that E2 > 2*E1^2
            % and the constraint that E3 > 3*E2^2/(2*E1)
            c(end+1)=(2+EPSTOL)*E1j(1)^2-E2j(1);
            c(end+1)=(3/2+EPSTOL)*E2j(1)^2/E1j(1)-E3j(1);
        end
        temp = prod(E3j);
        temp = temp/(factorial(3))^(NumMAPs-1);
        c(end+1) = temp/E(3) -2;
        c(end+1) = 0.5 - temp/E(3);
        for j = 2:NumMAPs
            c(end+1) =  1/3*E1j(j)/(E2j(j)-2*E1j(j)^2)*E3j(j)-1/2*E2j(j)^2/(E2j(j)-2*E1j(j)^2) - 1/(1e+16);
            c(end+1) =  1e-16 - 1/3*E1j(j)/(E2j(j)-2*E1j(j)^2)*E3j(j)-1/2*E2j(j)^2/(E2j(j)-2*E1j(j)^2);
            
        end
        
    end
    function f=objfun(x)
        [E1j,E3j]=xtopar(x);
        BCj=ones(1,size(BCLags,1));
        
        [compMAP, subMaps, err] = kpcfit_sub_compose(E1j, SCVj, E3j, G2j);
        if err ~= 0
            f=max(2*fold, 10^6);
            return
        end
        
        % scale MAP to match E(1)
        compMAP=map_scale(compMAP, E(1));
        
        for indexL=1:size(BCLags,1)
            BCj(indexL)=(map_joint(compMAP,BCLags(indexL,:),[1,1,1]));
        end
        
        f=norm(BC-BCj,2)/NBC;
        
        if isnan(f)
            f=2*fold;
        else
            fold=f;
        end
    end

end
