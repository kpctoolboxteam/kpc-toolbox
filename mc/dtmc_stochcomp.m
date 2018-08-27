function [S,P11,P12,P21,P22]=dtmc_stochcomp(P,I)
% function [S,P11,P12,P21,P22]=dtmc_stochcomp(P,I) - stochastic complementation in DTMCs
if nargin==1 
    I=1:ceil(length(P)/2);
end
Iall = 1:length(P);
Ic = Iall(~ismember(Iall,I)); % slightly faster than setdiff
P11=P(I,I);
P12=P(I,Ic);
P21=P(Ic,I);
P22=P(Ic,Ic);
S2 = eye(size(P22))-P22;
S=P11+P12*(S2 \ P21);
end
