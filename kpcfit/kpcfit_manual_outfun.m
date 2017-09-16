    function stop = kpcfit_manual_outfun(x,optimValues,state)
        stop = false;
        if ~isnan(optimValues.fval) & ~isnan(x)
            lgkx = x;
        else
            stop = true;
            %            fprintf('optimization stop: numerical difficulties detected\n')
        end
        if stagnval == 0
            stagnval = optimValues.fval;
        end
        delta = abs(optimValues.fval-stagnval)/stagnval;
        if delta < 0.01
            stagniter = stagniter + 1;
            if stagniter == 100
                %                fprintf('optimization stop: stagnation detected\n')
                stop = true;
            end
        end
        stagnval = optimValues.fval;
    end