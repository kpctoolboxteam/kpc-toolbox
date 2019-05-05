function [alpha,T,PHR]=map2ph(MAPIN)
% [MAPOUT,A,pi]=map2ph(MAPIN) - Returns a PH distribution
%
%  Input:
%  MAPIN: a MAP in the form of {D0,D1}
%
%  Output:
%  alpha: entry probability vector of the PH distribution
%  T: subgenerator of the PH distribution
%  PHR: PH-renewal process with distribution (pi,A) in (D0,D1) notation

PHR=MAPIN;
T=MAPIN{1};
alpha=dtmc_solve(-inv(MAPIN{1})*MAPIN{2});
PHR{2}=MAPIN{2}*ones((length(MAPIN{2})),1)*pi;
end

