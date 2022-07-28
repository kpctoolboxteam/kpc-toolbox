function scv=trace_scv(S)
% [scv]=trace_scv(S)
%
% DESCRIPTION
% Compute the squared coefficient of variation for trace S

scv = var(S)/mean(S)^2;
end