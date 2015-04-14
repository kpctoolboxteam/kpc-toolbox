function [BiCov,BiCovLags]=trace_bicov(S,GRID)
% compute the bicovariance of the trace
BiCovLags=[];
for i=GRID
    for j=GRID
        BiCovLags(end+1,:)=[1,i,j];
    end
end
BiCov = [];
for i=1:size(BiCovLags,1)
    BiCov(end+1)=trace_joint(S,BiCovLags(i,:),[1,1,1]);
end

end  
