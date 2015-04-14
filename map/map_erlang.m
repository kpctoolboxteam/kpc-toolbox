function MAP=map_erlang(MEAN,k)
% MAP=map_erlang(MEAN,k) - Fit an Erlang-k process as a MAP
%
%  Input: 
%  MEAN: mean inter-arrival time of the process
%  k: number of phases
%       
%  Output: 
%  MAP: a MAP in the form of {D0,D1}
%
%  Examples:
%  - MAP=map_erlang(2,3) an Erlang-3 process with mean 2
%

mu=k/MEAN;

MAP={zeros(k),zeros(k)};
% compute D0
for i=1:k-1
    MAP{1}(i,i+1)=mu;
end
% compute D1
MAP{2}(k,1)=mu;
% normalize D0 diagonal
MAP=map_normalize(MAP);
end