function [PH] = kpcfit_ph_auto(E,options)
% [PH] = kpcfit_ph_auto(E, options)
%
% DESCRIPTION
% Automatic fit a phase-type (PH) distribution using exact or approximate
% moment matching. Approximations are based on the Kronecker Product
% Composition method (KPC).
%
% INPUT
% E        - vector of moments, E(k) = E[X^k]
% options  - ph fitting
%
% OUTPUT
% PH{i,1}: the i-th fitted PH described as a PH-renewal process specified
%          in the (D0,D1) notation. The (alpha,T) representa
%
% EXAMPLE
%
% S; % trace
% for k=1:16 E(k) = mean(S.^k); end
% options = kpcfit_ph_options(E,'MinNumStates',2,'MaxNumStates',8);
% [PH]=kpcfit_ph_auto(E,options)
% kpcfit_ph_summary
% i=1; [alpha,T]=map2ph(PH{i,1}) % (alpha,T) distribution of i-th fitted PH

if options.Verbose
    fprintf('kpcfit_ph: version %s\n',kpcfit_version())
    fprintf('kpcfit_ph: type "help kpcfit_ph_options" for options\n')
end

%% initialization
% scale moments such that E[X]=1
Eunscaled = E;
Escale = [];
for k = 1:length(E)
    Escale(k) = E(1).^k;
end
E = (E(:) ./ Escale(:))';

%% determine best exact fit -- Prony's method
if options.Verbose
    fprintf('\nkpcfit_ph: starting exact fitting methods\n')
end
PH_EXACT = kpcfit_ph_exact(E, options);

%% determine best approximate fit -- moment space
if options.Verbose
    fprintf('\nkpcfit_ph: starting approximate fitting method -- moment space (KPC-based)\n')
    fprintf('\t\t\tRUN\t\tORD\tDIST\t\tRESULT')
end
PH_APX_MS = {};

tic;
for J = max([2,ceil(log2(options.MinNumStates))]):ceil(log2(options.MaxNumStates)) % for all number of MAPs to compose
    % initialize
    mindist = Inf;
    best = 0;
    bestidx = 0;
    
    % random runs
    randomruns = ceil(2*options.Runs/3); % number of runs with random initial point
    for run = 1:options.Runs
        if options.Verbose
            if run > randomruns
                fprintf(1,'\n\t\t\t%dr\t',run); %run number
            else
                fprintf(1,'\n\t\t\t%d\t',run); %run number
            end
            fprintf(1,'\t%d',2^J); % number of states
        end
        %try
        if run > randomruns && best > 0
            x0 = PH_APX_MS{bestidx,3}; %x0 = x0 * mean(x0)*1e-2*rand(size(x0));
        else
            x0 = []; % automatic
        end
        [APX, ~, x] = kpcfit_ph_search(E,J,options,x0);
        
        if map_isfeasible(APX)
            PH_APX_MS = store_result(PH_APX_MS,APX);
        else
            if options.Verbose
                fprintf(1,'\t\t\t\tx infeasible result.');
            end
        end
        %catch
        %    fprintf(1,'\t\t\t\tx optimization solver failed.');
        %end
    end
end

%% determine best approximate fit -- parameter space

SCV = (E(2)-E(1)^2)/E(1)^2;
if SCV >= 1
    if options.Verbose
        fprintf('\n\nkpcfit_ph: starting approximate fitting method -- hyper-exponential parameter space (non-KPC-based)\n')
    end
    PH_APX_PS = {};
    if options.Verbose
        fprintf('\t\t\tRUN\t\tORD\tDIST\t\tRESULT')
    end
    for J = max([2,ceil(log2(options.MinNumStates))]):ceil(log2(options.MaxNumStates)) % for all number of MAPs to compose
        % initialize
        mindist = Inf;
        best = 0;
        bestidx = 0;
        
        % random runs
        randomruns = ceil(2*options.Runs/3); % number of runs with random initial point
        for run = 1:options.Runs
            if options.Verbose
                fprintf(1,'\n\t\t\t%d\t',run); %run number
                fprintf(1,'\t%d',2^J); % number of states
            end
            try
                [APX,~,x] = psfit_hyperexp_weighted(E,2^J);
                
                if map_isfeasible(APX)
                    PH_APX_PS = store_result(PH_APX_PS,APX);
                else
                    if options.Verbose
                        fprintf(1,'\t\t\t\tinfeasible result.');
                    end
                end
            catch
                if options.Verbose
                    fprintf(1,'\t\t\t\tx optimization solver failed.');
                end
            end
        end
    end
else
    fprintf('\n\nkpcfit_ph: starting approximate fitting method -- aph parameter space (non-KPC-based)\n')
    PH_APX_PS = {};
    fprintf('\t\t\tRUN\t\tORD\tDIST\t\tRESULT')
    for J = max([2,ceil(log2(options.MinNumStates))]):ceil(log2(options.MaxNumStates)) % for all number of MAPs to compose
        % initialize
        mindist = Inf;
        best = 0;
        bestidx = 0;
        
        % random runs
        randomruns = ceil(2*options.Runs/3); % number of runs with random initial point
        for run = 1:options.Runs
            fprintf(1,'\n\t\t\t%d\t',run); %run number
            fprintf(1,'\t%d',2^J); % number of states
            try
                [APX,~,x] = psfit_aph_weighted(E,2^J);
                
                if map_isfeasible(APX)
                    PH_APX_PS = store_result(PH_APX_PS,APX);
                else
                    fprintf(1,'\t\t\t\tinfeasible result.');
                end
            catch
                if options.Verbose
                    fprintf(1,'\t\t\t\tx optimization solver failed.');
                end
            end
        end
    end
end

%% ensure that all returned results match at least the mean
nExactResults = size(PH_EXACT,1);
nApproxMS = size(PH_APX_MS,1);
%nApproxPS = size(PH_APX_PS,1);

for j = 1:size(PH_EXACT,1)
    PH{j,1} = map_scale(PH_EXACT{j},Eunscaled(1));
    PH{j,2} = dist_fun(E,map_moment(PH_EXACT{j},1:length(E)));
    PH{j,3} = [];
    PH{j,4} = 'exact';
end
for j = 1:size(PH_APX_MS,1)
    PH{nExactResults + j,1} = map_scale(PH_APX_MS{j,1},Eunscaled(1));
    PH{nExactResults + j,2} = dist_fun(E,map_moment(PH_APX_MS{j,1},1:length(E)));
    PH{nExactResults + j,3} = PH_APX_MS{j,3};
    PH{nExactResults + j,4} = 'approx_moment_space';
end
for j = 1:size(PH_APX_PS,1)
    PH{nExactResults + nApproxMS + j,1} = map_scale(PH_APX_PS{j,1},Eunscaled(1));
    PH{nExactResults + nApproxMS + j,2} = dist_fun(E,map_moment(PH_APX_PS{j,1},1:length(E)));
    PH{nExactResults + nApproxMS + j,3} = PH_APX_PS{j,3};
    PH{nExactResults + nApproxMS + j,4} = 'approx_param_space';
end

if options.Verbose
    fprintf(sprintf('\n\nReturned %d PH distribution(s) by exact moment matching.',nExactResults));
    fprintf(sprintf('\nReturned %d PH distribution(s) by approximate moment matching (moment space).',size(PH_APX_MS,1)));
    fprintf(sprintf('\nReturned %d PH distribution(s) by approximate moment matching (parameter space).',size(PH_APX_PS,1)));
    fprintf('\n\n');
end
return

    function dist = dist_fun(Eref,Eapx)
        % weight function
        w = (log10(Eref(end)).^(1/length(Eref))).^-(1:length(Eref));
        dist = w*abs( log10(Eref) - log10(Eapx) )'; %%does not work!
    end

    function PH_APX_MS = store_result(PH_APX_MS,APX)
        dist = dist_fun(E, map_moment(APX,1:length(E)));
        if isnan(dist)
            if options.Verbose
                fprintf(1,'\t-------- \tx no solution found.',dist);
            end
        else
            if options.Verbose
                fprintf(1,'\t%6.6f',dist);
            end
            if min([Inf,dist]) < mindist
                %            fprintf(1,'+',length(PH_APX_MS));
                if best == 0
                    PH_APX_MS{end+1,1} = APX;
                    if length(APX{1}) > 2^J
                        if options.Verbose
                            fprintf(1,sprintf('\to PH(%d) saved.  \t elapsed: %0.3fs [order automatically increased to match second moment]',length(APX{1}),toc));
                        end
                    else
                        if options.Verbose
                            fprintf(1,sprintf('\to PH(%d) saved.  \t elapsed: %0.3fs',length(APX{1}),toc));
                        end
                    end
                else
                    PH_APX_MS{end,1} = APX;
                    if options.Verbose
                        fprintf(1,sprintf('\t+ PH(%d) updated. \t elapsed: %0.3fs',length(APX{1}),toc));
                    end
                end
                best = run;
                bestidx = size(PH_APX_MS,1);
                mindist = min([Inf,dist]);
                PH_APX_MS{end,2} = dist;
                PH_APX_MS{end,3} = x;
            else
                if options.Verbose
                    fprintf(1,sprintf('\t\t\t\t\t\t elapsed: %0.3fs',toc));
                end
            end
        end
    end
end


