function pos=mapqn_hashstate(MAPQN,NVEC,KVEC)
% POS=mapqn_hashstate(MAPQN,NVEC,KVEC) - Determine position of a state
% (NVEC,KVEC) in the equilibrium probability vector of the model
%
%  Input: 
%  MAPQN: data structure summarizing the MAP network (see mapqn_ezsolve)
%  NVEC: distribution of jobs (see mapqn_ezsolve)
%  KVEC: distribution of active phases (see mapqn_ezsolve)
%       
%  Output: 
%  POS: position of state (NVEC, KVEC)
%
%  Examples:
%  - [XN,QN,UN,p,MAPQN]=mapqn_ezsolve({map_erlang(1,2);map_erlang(1,2)},4), 
%    p(mapqn_hashstate(MAPQN,[3,1],[2,2])) is the probability of the state
%    where queue 1 has 3 enqueued jobs and both erlang processes are in 
%    phase 2
%
% MAP Queueing Networks Toolbox
% Version 1.0 	 15-Apr-2008
SS=MAPQN.SS;
KK=MAPQN.KK;
sspos=matchrow(SS,NVEC);
kkpos=matchrow(KK,KVEC);
pos=(sspos-1)*size(KK,1)+kkpos;
    function feasiblerows = matchrow(Matrix, row)
        feasiblerows=1:size(Matrix,1);
        for col=1:length(row)
            A=Matrix(feasiblerows,col)==row(col);
            t=find(A);
            feasiblerows=feasiblerows(t);
            if length(feasiblerows)==1 && sum(Matrix(feasiblerows,:)==row)==size(Matrix,2)
                return
            end
        end
        feasiblerows=-1;
    end
end
