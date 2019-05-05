function stop = psfit_aph_weighted_outfun(x, optimValues, state)
global MAXTIME;
global MAXCHECKITER;
global T0;

stop = false;
if strcmpi(state,'iter')
    if mod(optimValues.iteration,MAXCHECKITER)==0 && optimValues.iteration>1
        reply = input('Do you want more? Y/N [Y]: ', 's');
        if isempty(reply)
            reply = 'Y';
        end
        if strcmpi(reply,'N')
            stop=true;
        end
    end
    if toc(T0)>MAXTIME
        fprintf('Time limit reached. Aborting.\n');
        stop = true;
    end
end
end