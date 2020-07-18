function [bestMAP,fac,fbc,subMAPs,otherMAPs, otherFACs, otherFBCs, otherSubMAPs]=kpcfit_manual(NumMAPs,E,AC,ACLags,BC,BCLags,varargin)
%
% DESCRIPTION
% Fitting of trace S into a Markovian Arrival Process based on Kronecker
% Product Composition (KPC) theory using manual parameter specification.
%
% OPTION LIST
% 'OnlyAC'    - boolean, if true fitting is performed on moments and
% 'MaxIterAC'  - integer, maximum number of iterations for a single autocorrelation fitting run
% 'MaxIterBC'  - integer, maximum number of iterations for a single bicovariance fitting run
% 'MaxRunsAC'  - maximum number of autocorrelation fitting runs
% 'MaxRunsBC'  - maximum number of bicovariance fitting runs
% 'MaxResAC'   - maximum number of autocorrelation fitting further analyzed for bicovariance fitting
warning off
%% options
OptionNames = [
    'MaxIterAC ';
    'MaxIterBC ';
    'MaxRunsAC ';
    'MaxRunsBC ';
    'MaxResAC  ';
    'OnlyAC    ';
    'AnimateAC ';
    'MaxRetMAPs';
    'ParallelBC';
    ];

OptionTypes = [
    'numeric';
    'numeric';
    'numeric';
    'numeric';
    'numeric';
    'numeric';
    'numeric';
    'numeric';
    'numeric'];

OptionValues = [];
for i = 1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end
% Default settings
options.OnlyAC    = false; % by default fit also bicovariances
options.AnimateAC = false; % by default fit also bicovariances
options.MaxIterAC = 300; % iterate up to 300 times to fit autocorrelations
options.MaxIterBC =  10; % iterate up to 10 times to fit bicorrelations
options.MaxRunsAC =  50; % maximum number of runs for AC fitting
options.MaxRunsBC =  30; % maximum number of runs for AC fitting
options.MaxResAC  =  min([options.MaxRunsAC,10]); % maximum number of values returned for AC fitting
options.MaxRetMAPs = 1;  % maximum number of MAPs returned
options.ParallelBC = 0; % parallelize BC runs


% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);


%% fit autocorrelations
disp(' ');
[resSCV,resG2,fobjAC] = kpcfit_sub_acfit(E,AC,ACLags,NumMAPs,options.MaxIterAC,options.MaxRunsAC,options.MaxResAC,options.AnimateAC);

%% fit bicovariances
disp(' ');
resE1 = cell(length(resSCV),1);
resE3 = cell(length(resSCV),1);
fobjBC = zeros(length(resSCV),1);

if options.OnlyAC == false    
    if options.ParallelBC > 0
        parfor i = 1:length(resSCV)
            fprintf(1,'fitting: processing acfit subresult %d of top %d\n',i,length(resSCV));
            [E1j,E3j,foBC] = kpcfit_sub_bcfit(E,resSCV{i},resG2{i},BC,BCLags,options.MaxIterBC,options.MaxRunsBC);
            resE1{i}  = E1j;
            resE3{i}  = E3j;
            fobjBC(i) = foBC;
        end
    else
         for i = 1:length(resSCV)
            fprintf(1,'fitting: processing acfit subresult %d of top %d\n',i,length(resSCV));
            [E1j,E3j,foBC] = kpcfit_sub_bcfit(E,resSCV{i},resG2{i},BC,BCLags,options.MaxIterBC,options.MaxRunsBC);
            resE1{i}  = E1j;
            resE3{i}  = E3j;
            fobjBC(i) = foBC;
        end
    end
elseif options.OnlyAC == true
    for i = 1:length(resSCV)
        resE1{i}  = ones(1,NumMAPs);
        resE3{i}  = (3/2+0.01).*(1+resSCV{i}).^2;
        fobjBC(i) = -1;
    end
else
    error('value of OnlyAC must be true or false')
end

%% determine best result
[v,ind] = sort(fobjBC,1,'ascend');

% truncate bicorrelations when it gets negative
tag=[];
for index=1:size(BCLags,1)
    if find(diff(BCLags(index,:))<0)
        tag(end+1,1)=index;
    end
end
BCLags(tag,:)=[];
BC(tag)=[];

%% compose intermediate results into final MAPs
composedMAPs = 0;
MAPs = cell(1,0);
subs = cell(1,0);
FACs = [];
FBCs = [];

for k=1:length(ind)
    index = ind(k);

    COMPOSE_VERBOSE=1;
    [newMAP,newSubMAPs,errorCode]=kpcfit_sub_compose(resE1{index},resSCV{index},resE3{index},resG2{index}, COMPOSE_VERBOSE);
    if errorCode ~= 0 || isempty(newMAP)
        fprintf("[%d] Discarded (errorCode=%d).\n", k, errorCode);
        continue
    end
    
    % scale MAP to match mean exactly
    newMAP=map_scale(newMAP,E(1));
    
    % Compute value of obj. functions for the resulting MAP
    [newfAC,newfBC]=evaluate_obj_function(newMAP);
    
    composedMAPs = composedMAPs + 1;

    MAPs{composedMAPs} = newMAP;
    FACs(composedMAPs) = newfAC;
    FBCs(composedMAPs) = newfBC;
    subs{composedMAPs} = newSubMAPs;
    
    if COMPOSE_VERBOSE > 0
        fprintf(1, "[%d] OK - fac= %f; fbc= %f\n", k, newfAC, newfBC);
    end

    if composedMAPs == options.MaxRetMAPs
        break
    end
end


%%
if composedMAPs == 0
    fprintf(2, "+++ KPC FAILED +++\n");
else
    bestMAP=MAPs{1};
    subMAPs=subs{1};
    fbc=FBCs(1);
    fac=FACs(1);
    fprintf(1, "Returned %d MAPs:\n", composedMAPs);
    fprintf(1, '1) fAC=%f, fBC=%f, SCV=%f, ACF(1)=%f, SKEW=%f\n', ...
        fac, fbc, map_scv(bestMAP), map_acf(bestMAP,1), map_skew(bestMAP));
    
    otherMAPs=cell(1,0);
    otherFBCs=[];
    otherFACs=[];
    otherSubMAPs=cell(1,0);
    
    
    for i=1:composedMAPs-1
        otherMAPs{end+1} = MAPs{1+i};
        otherFACs(end+1) = FACs(1+i);
        otherFBCs(end+1) = FBCs(1+i);
        otherSubMAPs{end+1} = subs{1+i};
        fprintf(1, '%d) fAC=%f, fBC=%f, SCV=%f, ACF(1)=%f, SKEW=%f\n', ...
            i+1, otherFACs(i), otherFBCs(i), map_scv(otherMAPs{i}), ...
            map_acf(otherMAPs{i},1), map_skew(otherMAPs{i}));
    end
    
    fprintf(1, '\n');
end



    function [objAC,objBC] = evaluate_obj_function (map)
        tSCV=(E(2)-E(1)^2)/E(1)^2;
        objAC=norm((AC-map_acf(map, ACLags)'),1)/norm(AC,2) + (map_scv(map) - tSCV)^2/tSCV^2;
        
        mapBC=ones(1,size(BCLags,1));
        for j=1:size(BCLags,1)
            mapBC(j)=map_joint(map,BCLags(j,:),[1,1,1]);
        end
        objBC = norm(BC-mapBC,2)/norm(BC,2);
    end


warning on
end
