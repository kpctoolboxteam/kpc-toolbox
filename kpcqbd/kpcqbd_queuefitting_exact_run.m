load('data/simulation_probs_avgtable_buildserver_1e10.mat')


clear PH;
clear scores;
clear PH_L;
clear exactfittingtoc;

repetitions = 3;
probInd = 1;
noOfPointsForFitting = 20; 
rho = 0.25;

K = [2,2];
J = [3,3];


scores = 10e3*ones(length(K),repetitions);
ARV = map_hyperexp(1/rho,5);

for k=1:length(K)
    for i=1:repetitions
        tic
        [PHcurr,scoreCurr,eflagcurr,xcurr,PH_Lcurr,x0curr] = kpcqbd_fit_exact(ARV,J(k),K(k),prob{probInd}(1:noOfPointsForFitting));
        exactfittingtoc(k,i)=toc;
        scores(k,i) = scoreCurr;
        currMinVec = min(scores,[],2);
        if  scores(k,i) <= currMinVec(k)
            for j=1:length(PHcurr)
                PH{k,j} = PHcurr{1,j};
            end
            x{k} = xcurr;
            x0{k} = x0curr;
            PH_L{k} = PH_Lcurr;
        end  
    end
end

scvs = [];
for i=1:length(PH_L)
    for j=1:length(PH_L{i})
        scvs(i,j) = map_scv(PH_L{i}{j});
    end
end

