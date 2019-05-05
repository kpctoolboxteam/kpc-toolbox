function PH_EXACT = kpcfit_ph_exact(E,options)
% Not for direct call, auxiliary function.
PH_EXACT = {};

SCV = ( E(2) - E(1)^2 ) / E(1)^2;
if SCV > 1
    fprintf(sprintf('kpcfit_ph: HIGHER variability than an exponential (var/mean^2 = %f)\n\n',SCV));
    fprintf('kpcfit_ph: starting exact hyper-exponential fitting method (Prony\''s method)\n');
    for n = 2:options.MaxNumStates
        if length(E) < 2*n-1
            fprintf(sprintf('kpcfit_ph: not enough moments given in input to fit hyper-exp(%d)\n',n))
            break
        end
        PH = kpcfit_ph_prony(E, n);
        if map_isfeasible(PH) > 0
            PH_EXACT{end+1} = PH;
            fprintf(sprintf('\t\t\thyper-exp(%d): feasible, matched exactly %d moments. result saved.\n',n,2*n-1))
        else
            fprintf(sprintf('\t\t\thyper-exp(%d): infeasible to fit exactly.\n',n))
            break
        end
    end
elseif SCV < 1
    fprintf(sprintf('kpcfit_ph: LOWER variability than an exponential (var/mean^2 = %f)\n\n',SCV))
    n = 1;
    while 1/n > SCV
        n = n + 1;
    end
    fprintf(sprintf('kpcfit_ph: exact fitting of E[X^2] requires at least %d states\n',n))
    if options.MinExactMom >= 2 & options.MaxNumStates < n
        fprintf(sprintf('kpcfit_ph: impossible to fit exactly E[X^2] with MaxNumStates = %d, increasing to MaxNumStates = %d.\n',options.MaxNumStates,n))
        return
    end
    
    if n == 2
        fprintf(sprintf('kpcfit_ph: attempting PH(2) fitting method\n',n));
        PH = map2_fit(E(1),E(2),E(3),0);
        if isempty(PH)
            PH = map2_fit(E(1),E(2),-1,0);
            if map_isfeasible(PH)
                PH_EXACT{end+1} = PH;
                fprintf(sprintf('\t\t\tph(2): feasible, matched exactly 2 moments. result saved.\n',n))
            else
                error('anomalous set of moments, please check.');
            end
        else
            PH_EXACT{end+1} = PH;
            fprintf(sprintf('\t\t\tph(2): feasible, matched exactly 3 moments. result saved.\n',n))
        end
    elseif (SCV > 1/n - kpcfit_tol) && (SCV < 1/n + kpcfit_tol)
        ERL = map_erlang(E(1),n); % exponential PH-renewal
        if norm(E - map_moment(ERL,1:length(E))) < kpcfit_tol
            fprintf(sprintf('kpcfit_ph: erlang moment set. fitted erlang-%d. result saved.\n',n))
            PH_EXACT{end+1} = ERL;
        end
    else
        %        fprintf('kpcfit_ph: no exact fitting method available for this moment set.\n');
        maxorder = options.MaxNumStates;
        fprintf(sprintf('kpcfit_ph: fitting APH distribution (best effort, max order = %d).\n',maxorder));
            warning off
        PH = aph_fit(E(1),E(2),E(3),maxorder);
        if map_isfeasible(PH)
%            fprintf(sprintf('kpcfit_ph: fitted APH(%d) distribution.\n',length(PH{1})));
            aph_matched = 0;
            for k=1:(2*length(PH{1})-1)
                if abs(E(k)-map_moment(PH,k)) < kpcfit_tol*map_moment(PH,k)
                    aph_matched = aph_matched + 1;
                end
            end
            fprintf(sprintf('\t\t\t      aph(%d): feasible, matched exactly %d moments. result saved.\n',length(PH{1}),aph_matched))
            PH_EXACT{end+1} = PH;
        else
            fprintf('kpcfit_ph: cannot fit APH distribution.\n');
        end
    end
    
else
    fprintf(sprintf('kpcfit_ph: SAME variability than an exponential (var/mean^2 = %f)\n\n',SCV))
    EXP = {[-1/E(1)],[1/E(1)]}; % exponential PH-renewal
    if norm(E - map_moment(EXP,1:length(E))) < kpcfit_tol
        fprintf('kpcfit_ph: exponential moment set. fitted exponential. result saved.\n')
    end
    PH_EXACT{end+1} = EXP;
end

end