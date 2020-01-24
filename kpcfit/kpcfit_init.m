function trace = kpcfit_init(S, varargin)
% trace = kpcfit_init(S,'option1','val1','option2','val2',...)
%
% DESCRIPTION
% Prepare trace S for KPC fitting
%
% OUTPUT
% trace.S - original trace
% trace.E - moments of the trace, e.g., E(3) is E[X^3]
% trace.ACFull - autocorrelation coefficients of the trace. Each coefficient is computed at least from 10 samples.
% trace.AC - autocorrelation coefficients selected for fitting
% trace.BC - bicovariance coefficients selected for fitting (BCFull is too expensive to compute)
% trace.BCGridLags - bicovariance lags (2-dimensional sampling square grid)
% trace.BCLags - bicovariance lags (1-dimensional side of the sampling square grid)
%
% OPTION LIST
% 'ACLags'     - sets lags of of the autocorrelation coefficients
% 'BCGridLags' - sets lags of of the sampling square grid
% 'Smooth' - smooths autocorrelation coefficients using the specified window size

%% options
OptionNames = [
    'ACLags    ';
    'BCGridLags';
    'Smooth    '
    ];

OptionTypes = [
    'numeric';
    'numeric';
    'numeric'];

OptionValues = [];
for i = 1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end

% Default settings
n = length(S);
options.MaxMoment = 3; % fit 3 moments (this value is currently unused)
nMinSupportAC = 10; % mininum number of points needed to estimate an AC coefficient
options.ACLags = unique(logspacei(1,ceil(n/nMinSupportAC),500)); % AC lags used for fitting
options.BCGridLags = unique(logspacei(1,max(options.ACLags),5)); % lags used to generate square BC grid for fitting
options.Smooth = 0;
options.LogSmooth = 0;
% Parse optional parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

%% trace preprocessing
fprintf(1,'init: computing moments from the trace\n');
% moments
E=zeros(1,options.MaxMoment);
for j=1:options.MaxMoment
    E(j)=mean(S.^j);
end

if options.Smooth > 0
    fprintf(1,'init: computing smoothed autocorrelations from the trace\n');
    if exist('smooth')==0
     warning('MATLAB:kpcfit_init:cftoolbox','curve fitting toolbox not installed, smooth function unavailable\n');
     fprintf(1,'init: computing autocorrelations from the trace\n');  
     AC=trace_acf(S,options.ACLags);
     ACFull=trace_acf(S,1:ceil(length(S)/nMinSupportAC));
    else
     AC=smooth(trace_acf(S,options.ACLags),options.Smooth);
     ACFull=smooth(trace_acf(S,1:ceil(length(S)/nMinSupportAC)),options.Smooth);
    end
else
    fprintf(1,'init: computing autocorrelations from the trace\n');
    AC=trace_acf(S,options.ACLags);
    ACFull=trace_acf(S,1:ceil(length(S)/nMinSupportAC));
end
% autocorrelations

%

% estimate cut point of AC
posmax = options.ACLags(find(abs(AC) < 10^-6, 1 ));
%fprintf(1,'init: autocorrelations beyond lag %d ignored (abs. val. <10^-6)\n',posmax);
if ~isempty(posmax)
    todel = find(options.ACLags>posmax);
    options.ACLags(todel)=[]; % lags of AC used for fitting
    AC(todel)=[];
    todel = options.BCGridLags>posmax;
    options.BCGridLags(todel)=[]; % lags of AC used for fitting
end
fprintf(1, "Using %d AC Lags\n", length(options.ACLags));


fprintf(1,'init: computing bicovariances from the trace\n');
[BC,BCLags]=trace_bicov(S,options.BCGridLags); % bispectrum coefficients

%% generate trace structure
trace = struct('S',S(:),'E',E,'AC',AC(:),'ACFull',ACFull(:),'ACLags',options.ACLags,'BC',BC,'BCGridLags',options.BCGridLags,'BCLags',BCLags);
end
