function stop = kpcfit_sub_acfit_outfun(x,optimValues,state,MaxIterAC,MaxTimeAC,tstart,f_best)
global lgkx;
global stagnval;
global stagniter;
stop = false;
if strcmpi(state,'iter')
    if mod(optimValues.iteration,MaxIterAC)==0 && optimValues.iteration>1
        if ( optimValues.fval >f_best)
            stop = true;
        end
    end
    telapsed = toc(tstart);
    if ( telapsed>MaxTimeAC && optimValues.iteration > MaxIterAC)
        fprintf(1,'acfit: time limit reached in autocorrelation fitting\n');
        stop = true;
    end
end
if ~isnan(optimValues.fval) & sum(isnan(x)) == 0
    lgkx = x; % last good known solution
else
    stop = true;
    %fprintf(' optimization halted: numerical difficulties\n')
end
if stagnval == 0
    stagnval = optimValues.fval;
end
delta = abs(optimValues.fval-stagnval)/stagnval;
if delta < 0.01
    stagniter = stagniter + 1;
    if stagniter == 100
        %fprintf(' optimization halted: stagnation ')
        stop = true;
    end
else
    stagniter = 0;
end
stagnval = optimValues.fval;
end