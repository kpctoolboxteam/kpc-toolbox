function responseSample = GG1ResponseTime(arrTrace,serviceTrace,sampleSize)

    queue = javaObject('java.util.LinkedList');

    arrLength = length(arrTrace);
    serLength = length(serviceTrace);

    timeToNextA = arrTrace(randi(arrLength));
    timeToNextS = serviceTrace(randi(serLength));

    queueLength = 0;
    time = 0;
    responseSample = [];

    for i=1:sampleSize
        queueLength = queue.size();
        if(queueLength == 0)
            % empty queue, jump to next arr
            time = time + timeToNextA;
            queue.add(time);          

            timeToNextA = arrTrace(randi(arrLength));
        else
            if(timeToNextA < timeToNextS)
                % next event is arr
                time = time + timeToNextA;
                queue.add(time);

                timeToNextS = timeToNextS - timeToNextA;
                timeToNextA = arrTrace(randi(arrLength));
            else
                % next event is dep
                time = time + timeToNextS;
                arrivedAt = queue.remove();

                if i>sampleSize*0.1
                    responseSample(end+1) = time - arrivedAt;
                end

                timeToNextA = timeToNextA - timeToNextS;
                timeToNextS = serviceTrace(randi(serLength));
            end
        end
    end
end
