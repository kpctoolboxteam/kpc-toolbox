function stop = kpcfit_ph_outfun(x,optimValues,state)
global stagnval;
global stagniter;
global lgkx;
stop = false;
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
