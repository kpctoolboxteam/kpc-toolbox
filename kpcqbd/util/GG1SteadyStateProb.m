function prob = GG1SteadyStateProb(arrTrace,serviceTrace,sampleSize,maxState)

    arrLength = length(arrTrace);
    serLength = length(serviceTrace);

    timeToNextA = arrTrace(randi(arrLength));
    timeToNextS = serviceTrace(randi(serLength));

    queueLength = 0;

    prob = zeros(1,maxState);

    for i=1:sampleSize
       
        if(queueLength == 0)
            % empty queue, jump to next arr
            if i>1e6 && queueLength < maxState
                prob(queueLength+1) = prob(queueLength+1)+timeToNextA;
            end
            queueLength = queueLength +1;          

            timeToNextA = arrTrace(randi(arrLength));
        else
            if(timeToNextA < timeToNextS)
                % next event is arr
                if i>1e6 && queueLength < maxState
                    prob(queueLength+1) = prob(queueLength+1)+timeToNextA;
                end
                queueLength = queueLength + 1;

                timeToNextS = timeToNextS - timeToNextA;
                timeToNextA = arrTrace(randi(arrLength));
            else
                % next event is dep
                if i>1e6 && queueLength < maxState                 
                    prob(queueLength+1) = prob(queueLength+1)+timeToNextS;
                end
                queueLength = queueLength - 1;

                timeToNextA = timeToNextA - timeToNextS;
                timeToNextS = serviceTrace(randi(serLength));
            end
        end
    end
    
    prob = prob./sum(prob);

end
