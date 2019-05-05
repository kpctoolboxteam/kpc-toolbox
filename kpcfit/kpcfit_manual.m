function [MAP,fac,fbc,subMAPs]=kpcfit_manual(NumMAPs,E,AC,ACLags,BC,BCLags,varargin)
% [MAP,fac,fbc,subMAPs]=kpcfit_manual(NumMAPs,E,AC,ACLags,BC,BCLags,varargin)
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
    ];

OptionTypes = [
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
    for i = 1:length(resSCV)
        fprintf(1,'\nfitting: processing acfit subresult %d of top %d\n',i,length(resSCV));
        [E1j,E3j,foBC] = kpcfit_sub_bcfit(E,resSCV{i},resG2{i},BC,BCLags,options.MaxIterBC,options.MaxRunsBC);
        resE1{i}  = E1j;
        resE3{i}  = E3j;
        fobjBC(i) = foBC;
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
[v,ind] = sort(fobjBC,2,'ascend');
bestpos = ind(1);
fac = fobjAC(bestpos);
fbc = fobjBC(bestpos);

%% compose intermediate results into final MAP
[MAP,subMAPs]=kpcfit_sub_compose(resE1{bestpos},resSCV{bestpos},resE3{bestpos},resG2{bestpos});

%% scale MAP to match mean exactly
MAP=map_scale(MAP,E(1));
warning on
end