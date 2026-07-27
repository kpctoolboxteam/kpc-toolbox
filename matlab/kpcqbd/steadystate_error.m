function steadyStateScore = steadystate_error(ARV,PH_L,simProb)
    phord = PH_L;
    try
        pa_kpc = kpcqbd_solve(ARV,phord,length(simProb))';
    catch
        steadyStateScore = 1e10;
        return;
    end
    
    steadyStateScore = norm((pa_kpc(1:length(simProb))-simProb)./simProb);
end