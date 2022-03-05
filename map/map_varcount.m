function varc=map_varcount(MAP,tset)
n = length(MAP{1});
Q = map_infgen(MAP);
I = eye(n);
e = ones(n,1);
piq = map_pi(MAP);
lambda = 1/map_mean(MAP);
PRE = (lambda - 2*lambda^2 + 2*piq*MAP{2}*(ones(n,1)*piq-Q)^(-1)*MAP{2}*e);
for j=1:length(tset)
    POST = -2*piq*MAP{2}*(I-expm(Q*tset(j)))*(ones(n,1)*piq-Q)^(-2)*MAP{2}*e
    varc = PRE*tset(j) + POST;
end
end