function stop = kpcfit_sub_bcfit_outfun(x,optimValues,state,MaxIterBC,MaxTimeBC,tstart,f_best)
%numIterations = optimValues.iteration;
stop = false;
if strcmpi(state,'iter')
    if optimValues.iteration >= MaxIterBC && optimValues.iteration>1
        if ( optimValues.fval > f_best)
            stop = true;
        end
    end
    telapsed = toc(tstart);
    if (telapsed>MaxTimeBC && optimValues.iteration > MaxIterBC)
        fprintf('Time limit reached in moments fitting. Aborting.\n');
        stop = true;
    end
end
end