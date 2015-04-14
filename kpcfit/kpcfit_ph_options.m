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
% 'MinNumStates' - integer, minimum number of states - ** KPC returns only multiples of 2 **, default: 2
% 'MaxNumStates' - integer, maximum number of states - ** KPC returns only multiples of 2 **, default: 32
% 'MinExactMom'  - integer, minimum number of moments to be fitted exactly (best-effort), default: 2
% 'Runs'         - integer, maximum number of runs, default: 5
%
% EXAMPLE
% options = kpcfit_ph_options(E) % default options
% options = kpcfit_ph_options(E,'MinExactMom',2,'MaxNumStates',4,'Runs',1)

% Default options
options.('Runs') = 5;
options.('MinNumStates') = 2;
options.('MaxNumStates') = 32;
options.('MinExactMom') = 3;

ranges.('Runs') = [1,Inf];
ranges.('MinNumStates') = [2,Inf];
ranges.('MaxNumStates') = [2,Inf];
ranges.('MinExactMom') =  [1,length(E)];

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

if 2*options.('MaxNumStates')-1 > length(E)
    error('MATLAB:kpcfit_ph_options:insufficient_moments',sprintf('MaxNumStates of %d requires to supply at least %d moments.', options.('MaxNumStates'),2*options.('MaxNumStates')-1));
end

if options.('MinNumStates') > options.('MaxNumStates')
    warning('MATLAB:kpcfit_ph_options:max_lt_min','MaxNumStates < MinNumStates, fixed.');
    tmp = options.('MinNumStates');
    options.('MaxNumStates') = tmp;
    options.('MinNumStates') = options.('MaxNumStates');
end

end
