function acfc=map_acfc(MAP,kset,u)
% Autocorrelation of a MAP's counting process at lags in kset
% u = length of timeslot (ie timescale)

n = length(MAP{1});
Q = map_infgen(MAP);
I = eye(n);
piq = map_prob(MAP);
PRE = piq*MAP{2}*(I-expm(Q*u));
POST = (I-expm(Q*u))*(ones(n,1)*piq-Q)^(-2)*MAP{2}*ones(n,1);
vart=map_varcount(MAP,u);
for j=1:length(kset)
    acfc(j)=PRE * expm(Q*(kset(j)-1)*u) * POST/vart;
end
end
