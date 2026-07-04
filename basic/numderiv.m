function d = numderiv(f, x0, order)
% d = NUMDERIV(F,X0,ORDER) - ORDER-th derivative of a scalar function F,
% analytic in a neighbourhood of X0, evaluated at X0.
%
% Uses Cauchy's integral formula sampled on a small circle and evaluated by
% an FFT: for F analytic, F^(k)(x0) = k!/(2*pi*i) * oint F(z)/(z-x0)^(k+1) dz.
% Spectrally accurate for entire integrands; repo-native replacement for the
% numeric-derivative use of the third-party derivest routine.

m = 64;               % samples on the circle (>> order)
r = 0.5;              % radius
jj = (0:m-1)';
z  = x0 + r*exp(2i*pi*jj/m);
fz = zeros(m,1);
for k = 1:m
    fz(k) = f(z(k));
end
c = fft(fz)/m;        % c(k+1) ~ F^(k)(x0) * r^k / k!
d = real(c(order+1) * factorial(order) / r^order);
end
