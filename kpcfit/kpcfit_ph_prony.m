function PH=kpcfit_ph_prony(E,n)
f = factorial(0:(2*n-1));
m = kpcfit_hyper_charpoly(E,n);
theta=roots(m); % determine eigenvalues from characteristic poly
C = [];
for i= 1:n
    C(i,1:n) = f(i+1)*theta.^i;
end
M = (C\E(1:n)')'; % determine entry probabilities
PH{1}=diag(-1./theta);
PH{2}=-diag(-1./theta)*e(n)*M;
end