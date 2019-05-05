function [MAP,mu1,mu2,p]=map_hyperexp(MEAN,SCV,p)
% MAP=map_hyperexp(MEAN,SCV,p) - Fit a two-phase Hyperexponential process as a MAP
%
%  Input:
%  MEAN: mean inter-arrival time of the process
%  SCV: squared coefficient of variation of inter-arrival times
%  p: probability of being served in phase 1 (DEFAULT: p=0.99)
%
%  Output:
%  MAP: a MAP in the form of {D0,D1}
%
%  Examples:
%  - MAP=map_hyperexp(1,2) a two-phase Hyperexponential process where
%    phase 1 is selected with probability 0.99
%  - MAP=map_hyperexp(1,2,0.2) a two-phase Hyperexponential process where
%    phase 1 is selected with probability 0.2
%  - map_isfeasible(map_hyperexp(1,25,0.5)) a two-phase Hyperexponential
%    with SCV=25 and p=0.5 does not exist
%

if nargin<3
    p=0.99;
end
% there are two possible solutions, try with the first root
E2=(1+SCV)*MEAN^2;
Delta=-4*p*MEAN^2+4*p^2*MEAN^2+2*E2*p-2*E2*p^2;
mu2=(-2*MEAN+2*p*MEAN+sqrt(Delta))/(E2*p-2*MEAN^2);
mu1=mu2*p/(p-1+MEAN*mu2);
MAP={[-mu1,0;0,-mu2],[mu1*p,mu1*(1-p);mu2*p,mu2*(1-p)]};
if map_isfeasible(MAP) % if the solution is feasible
    return
else % if the solution is not feasible go for the second root
    mu2=(-2*MEAN+2*p*MEAN-sqrt(Delta))/(E2*p-2*MEAN^2);
    mu1=mu2*p/(p-1+MEAN*mu2);
    MAP={[-mu1,0;0,-mu2],[mu1*p,mu1*(1-p);mu2*p,mu2*(1-p)]};
    if p>1e-6 && ~map_isfeasible(MAP) % if the solution is not feasible try to decrease p        
        [MAP,mu1,mu2,p] = map_hyperexp(MEAN,SCV,p/10);
    else
        MAP = {};
        return
    end
end
end
