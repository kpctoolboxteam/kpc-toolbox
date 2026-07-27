function v = map_count_var(map,t)
% Computes the variance of the counting process, at resolution t, for the
% given MAP.
% Input:
% - map: Marovian Arrival Process (symbolic or numeric)
% - t: the period considered for each sample of the counting process
% Output:
% - v: varianze of arrivals in (0,t]
% Reference: [He and Neuts, Markov chains with marked transitions, 1998]
% Verified special case of MMPP(2) with [Andresen and Nielsen, 1998].

n = size(map{1},1);

D = map{1} + map{2};
if map_issym(map)
    I = sym(eye(n));
    e = sym(ones(n,1));
    v = sym(zeros(length(t),1));
else
    I = eye(n);
    e = ones(n,1);
    v = zeros(length(t),1);
end
theta = map_prob(map);
tmp = (e * theta - D)^(-1);

% arrival rate
l = theta * map{2} * e;
c = theta * map{2} * tmp;
d = tmp * map{2} * e;
ll = theta * map{2} * e;
for i=1:length(t)
    v(i) = (ll-2*l^2 + 2*c*map{2}*e)*t(i) - 2*c*(I-expm(D*t(i)))*d;
end

end