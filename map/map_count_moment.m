function [M] = map_count_moment(MAP, t, order)
% Computes power moments of counts, at resolution t, of a MAP.
% INPUT
% - MAP: Markovian Arrival Process
% - t: resolution
% - order: orders of the moments to compute
% OUTPUT
% - M: power moments of counts

n = size(MAP{1},1);
theta = map_prob(MAP);

if map_issym(MAP)
    e = sym(ones(n,1));
    M = sym(zeros(length(order),1));
else
    e = ones(n,1);
    M = zeros(length(order),1);
end

if map_issym(MAP) || max(order) > 4
    % symbolic derivative
    z = sym('z');
    MZ = theta*expm(MAP{1}*t+MAP{2}*exp(z)*t)*e;
    for i = 1:length(order)
        M(i) = subs(diff(MZ,z,order(i)),z,0);
    end
else
    % numerical derivative
    for i = 1:length(order)
        M(i) = derivest(@(z) mgfunc(z), 0, 'deriv', order(i));
    end
end

    function r = mgfunc(z)
        for j=1:length(z) %derivest passes a vector
            r(j) = theta*expm(MAP{1}*t+MAP{2}*exp(z(j))*t)*e;
        end
    end

end