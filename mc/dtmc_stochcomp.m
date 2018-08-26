function [S,P11,P12,P21,P22]=dtmc_stochcomp(P,I)
% function [S,P11,P12,P21,P22]=dtmc_stochcomp(P,I) - stochastic complementation in DTMCs
if nargin==1 
    I=1:ceil(length(P)/2);
end
%Ic=setdiff(1:length(P),I);
Iall = 1:length(P);
Ic = Iall(~ismember(Iall,I)); % slightly faster than setdiff
P11=P(I,I);
P12=P(I,Ic);
P21=P(Ic,I);
P22=P(Ic,Ic);
S2 = eye(size(P22))-P22;
%zero_row = find(sum(S2,2)==0);
%zero_col = find(sum(S2,1)==0);
%S2(zero_row,zero_row) = -eye(length(zero_row)); % can this be replaced by []?
%S2(zero_col,zero_col) = -eye(length(zero_col));
S=P11+P12*(S2 \ P21);
end
