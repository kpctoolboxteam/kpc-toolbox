function [fit] = mmpp2_fit_count_theoretical(map,t1,t2,tinf)
% Fits the theoretical characteristics of a MAP(n) with a MMPP(2).
if nargin<4
t1 = 1;
t2 = 10;
tinf = 1e8;
end

% joint-process charactersitics
a = map_count_mean(map,t1)/t1;
bt1 = map_count_var(map,t1)/(a*t1);
bt2 = map_count_var(map,t2)/(a*t2);
binf = map_count_var(map,tinf)/(a*tinf);
mt2 = map_count_moment(map,t2,1:3);
m3t2 = mt2(3) - 3*mt2(2)*mt2(1) + 2*mt2(1)^3;

fit = mmpp2_fit_count(a, bt1, bt2, binf, m3t2, t1, t2);