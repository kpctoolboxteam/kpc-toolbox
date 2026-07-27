function options = kpcfit_ph_options(E,varargin)
% options = kpcfit_ph_options(E)
%
% DESCRIPTION
% Specify options for kpcfit_ph(E,...)
%
% INPUT
% E              - vector of moments, E(k) = E[X^k]
%
% OPTIONS
% 'MinNumStates' - integer, minimum number of states - ** KPC accepts only multiples of 2 **, default: 2
% 'MaxNumStates' - integer, maximum number of states - ** KPC accepts only multiples of 2 **, default: 32
% 'MinExactMom'  - integer, minimum number of moments to be fitted exactly (best-effort), default: 2
% 'MaxExactMom'  - integer, maximum number of moments to be fitted exactly (best-effort), default: 2
% 'MinApproxMom' - integer, minimum number of moments to be fitted approximately (best-effort), default: 2
% 'MaxApproxMom' - integer, maximum number of moments to be fitted approximately (best-effort), default: 2
% 'MaxRuns'      - integer, maximum number of runs, default: 5
%
% EXAMPLE
% options = kpcfit_ph_options('MinExactMom',2,'MaxExactMom',4,'MaxRuns',1)

% Default options
options.('MinNumStates') = 2;
options.('MaxNumStates') = 32;
options.('MinApproxMom') = 2;
options.('MaxApproxMom') = length(E);
options.('MinExactMom') = 2;
options.('MaxExactMom') = 2;

ranges.('MaxRuns') = [1,Inf];
ranges.('MinNumStates') = [2,Inf];
ranges.('MaxNumStates') = [Inf,Inf];
ranges.('MinApproxMom') = [1,length(E)];
ranges.('MaxApproxMom') = [1,length(E)];
ranges.('MinExactMom') =  [1,length(E)];
ranges.('MaxExactMom') =  [1,length(E)];
ranges.('MaxRuns') = [1,Inf];

% Parse optional parameters
if mod(length(varargin),2) > 0 % if odd number
    error('MATLAB:kpcfit_ph_options:odd_number_of_parameters','Odd number of input parameters, please specify a value after each option.');
end

for i=1:2:length(varargin)
    options.(varargin{i}) = varargin{i+1};
end

% Sanity checks
for i=1:2:length(varargin)
    if options.(varargin{i}) < ranges.(varargin{i})(1)
        options.(varargin{i}) = ranges.(varargin{i})(1);
        warning('MATLAB:kpcfit_ph_options:out_of_range_parameter',sprintf('The minimum value for option %s is %d, now fixed.', varargin{i}, ranges.(varargin{i})(1)));
    elseif options.(varargin{i}) > ranges.(varargin{i})(2)
        options.(varargin{i}) = ranges.(varargin{i})(2);
        warning('MATLAB:kpcfit_ph_options:out_of_range_parameter',sprintf('The maximum value for option %s is %d, now fixed.', varargin{i}, ranges.(varargin{i})(2)));
    end
end

if mod(options.('MinNumStates'),2^round(log2(options.('MinNumStates')))) > 0 % is multiple of 2?
    options.('MinNumStates') = 2^ceil(log2(ceil(options.('MinNumStates'))));
    warning('MATLAB:kpcfit_ph_options:not_power_of_two',sprintf('MinNumStates not a multiple of 2, fixed to %d.', options.('MinNumStates')));
end

if mod(options.('MaxNumStates'),2^round(log2(options.('MaxNumStates')))) > 0 % is multiple of 2?
    options.('MaxNumStates') = 2^ceil(log2(ceil(options.('MaxNumStates'))));
    warning('MATLAB:kpcfit_ph_options:not_power_of_two',sprintf('MaxNumStates not a multiple of 2, fixed to %d.', options.('MaxNumStates')));
end

if options.('MaxNumStates') > options.('MinNumStates')
    warning('MATLAB:kpcfit_ph_options:max_lt_min','MaxNumStates < MinNumStates, fixed swapping values.');
    tmp = options.('MinNumStates');
    options.('MaxNumStates') = tmp;
    options.('MinNumStates') = options.('MaxNumStates');
end

if options.('MinApproxMom') > options.('MaxApproxMom')
    warning('MATLAB:kpcfit_ph_options:max_lt_min','MaxApproxMom < MinApproxMom, fixed swapping values.');
    tmp = options.('MinApproxMom');
    options.('MaxApproxMom') = tmp;
    options.('MinApproxMom') = options.('MaxApproxMom');
end

if options.('MinExactMom') > options.('MaxExactMom')
    warning('MATLAB:kpcfit_ph_options:max_lt_min','MaxExactMom < MinExactMom, fixed swapping values.');
    tmp = options.('MinExactMom');
    options.('MaxExactMom') = tmp;
    options.('MinExactMom') = options.('MaxExactMom');
end

end
