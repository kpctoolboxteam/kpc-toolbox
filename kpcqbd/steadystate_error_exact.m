function steadyStateScore = steadystate_error_exact(ARV,PH_L,simProb)
    PH = PH_L{1};

    for i=2:length(PH_L)
        PH = map_kpc(PH,PH_L{i});
    end

    try
        [~,~,~,pa_exact] = kpcqbd_mapmap1(ARV,PH);
    catch
        steadyStateScore = 1e10;
        return;
    end

    noOfLevels = min(size(pa_exact,1),length(simProb));
    
    steadyStateScore = norm((sum(pa_exact(1:noOfLevels,:),2)'-simProb(1:noOfLevels))./simProb(1:noOfLevels));
end