function [MAP,fac,fbc,kpcMAPs,otherMAPs, otherFACs, otherFBCs, otherSubMAPs]=kpcfit_auto(trace, varargin)
% [MAP,fac,fbc,kpcMAPs] = kpcfit_auto(T,'option1',val1,'option2',val2,...)
%
% DESCRIPTION
% Automatic fitting of trace T into a Markovian Arrival Process based on
% the Kronecker Product Composition (KPC) method.
%
% INPUT
% T        - data structure returned by kpcfit_init
% 
% OUTPUT
% MAP          - fitted MAP
% fac          - best-found objective function value for autocorrelation fitting (smaller = better)
% fbc          - best-found objective function value for bicorrelation fitting (smaller = better)
% kpcMAPs      - MAP(2)s to were composed by Kronecker products to define MAP
%
% EXAMPLE
%  T = kpcfit_init(S)
%  MAP = kpcfit_auto(T,'NumMAPs',3,'MaxIterAC',10)
%
% OPTION LIST
% 'OnlyAC'     - if 1 fitting is performed on moments and autocorrelations
% 'NumStates'  - integer, equals log2(number of states of the final MAP), default: BIC heuristic
% 'NumMAPs'    - integer, equals log2(number of states of the final MAP), default: BIC heuristic
% 'MaxRunsAC'  - maximum number of autocorrelation fitting runs (run = optimization program execution)
% 'MaxRunsBC'  - maximum number of bicovariance fitting runs
% 'MaxIterAC'  - integer, maximum number of iterations for a single autocorrelation fitting run
% 'MaxIterBC'  - integer, maximum number of iterations for a single bicovariance fitting run
% 'MaxResAC'   - maximum number of autocorrelation fitting runs further considered in bicovariance fitting
% 'ACLags'     - autocorrelation lags to be fitted (e.g., 1:10)
% 'BCGridLags' - lags defining bicovariance grid to be fitted  (e.g., 1:5 gives a 25 points grid 1:5x1:5)
% 'MaxRetMaps' - maximum number of MAPs returned
% 'ParallelBC' - enable parallelization of BC runs
%
% REFERENCES
% [1] G.Casale, E.Z.Zhang, E.Smirni. Trace Data Characterization and Fitting 
%     for Markov Modeling, Elsevier Performance Evaluation, 67(2):61-79, 
%     Feb 2010.
%
% [2] G.Casale, E.Z.Zhang, E.Smirni. KPC-Toolbox: Simple Yet Effective Trace
%     Fitting Using Markovian Arrival Processes. in Proc. of QEST 2008, 
%     83-92, St.Malo, France, IEEE Press, September 2008. 
%

%% options
OptionNames = [
    'OnlyAC    ';
    'NumStates ';
    'NumMAPs   ';
    'MaxIterAC ';
    'MaxIterBC ';
    'MaxRunsAC ';
    'AnimateAC ';
    'MaxRunsBC ';
    'MaxResAC  ';    
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
    'numeric';
    'numeric';
    'numeric'];

OptionValues = [];
for i = 1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end

% Default settings
options.OnlyAC    = false; 
options.AnimateAC = false; 
options.NumMAPs   = []; 
options.NumStates = []; 
options.MaxIterAC = 300; % iterate up to 300 times to fit autocorrelations
options.MaxIterBC =  10; % iterate up to 10 times to fit bicorrelations
options.MaxRunsAC =  50; % maximum number of runs for AC fitting
options.MaxRunsBC =  30; % maximum number of runs for AC fitting
options.MaxResAC  =  min([options.MaxRunsAC,10]); % maximum number of values returned for AC fitting
options.MaxRetMAPs = 1; % maximum number of MAPs returned
options.ParallelBC = 0; % parallelize BC runs

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);
if mod(options.NumStates,2) > 0 
   error('error: requested %d states, but kpc-toolbox can return only a number of states that is a power of 2');
end
if isempty(options.NumMAPs)    
    options.NumMAPs = ceil(log2(options.NumStates)); 
end

disp('** KPC fitting algorithm initialized **');
%% fitting algorithm run parameterization
% order
if isempty(options.NumMAPs) % if MAP order not given, use automatic BIC selection
    disp('init: performing order selection');
    try
        options.NumMAPs = kpcfit_sub_bic(trace.ACFull,[2 4 8 16 32 64 128]);
    catch
        fprintf(2, "BIC order selection failed. Using 3 MAPs\n");
        options.NumMAPs = 3;
    end
end
fprintf(1,'init: kpc-toolbox will search for a MAP with 2^%d=%d states\n',options.NumMAPs,2^options.NumMAPs);
%% fitting
disp('fitting: running KPC-based fitting script');    
[MAP,fac,fbc,kpcMAPs,otherMAPs, otherFACs, otherFBCs, otherSubMAPs]=kpcfit_manual(options.NumMAPs, ...
                            trace.E, ...
                            trace.AC, ...
                            trace.ACLags, ...
                            trace.BC, ...
                            trace.BCLags, ...
                            'MaxIterAC',options.MaxIterAC, ...
                            'MaxRunsAC',options.MaxRunsAC, ...
                            'MaxResAC',options.MaxResAC, ...
                            'MaxRunsBC',options.MaxRunsBC, ...
                            'MaxIterBC',options.MaxIterBC, ...
                            'OnlyAC',options.OnlyAC, ...
                            'AnimateAC',options.AnimateAC, ...
                            'MaxRetMAPs', options.MaxRetMAPs, ...
                            'ParallelBC', options.ParallelBC);

%% display final moments and autocorrelations
disp('** KPC fitting algorithm completed ** ');
disp(' ');
disp('                  Moments Comparison       ');
disp('          Original Trace          Fitted MAP');
format long e
for k=1:length(trace.E)
    disp([trace.E(k) map_moment(MAP,k)]); 
end
disp(' ');
disp('                             Autocorrelation Comparison       ');
disp('          Lag                     Original Trace          Fitted MAP');
format long e
for k=1:min([10,length(trace.ACLags)])
    disp([trace.ACLags(k) trace.AC(k) map_acf(MAP,k)]); 
end

end
