function m = hyper_charpoly(E,n)
E=[1,E(:)'];
f = factorial(0:(2*n-1));
for i = 1:n
    A(i,:) = E((n+i):-1:i)./f((n+i):-1:i);
end
A(n+1,:) = zeros(1,n+1); A(n+1,1) = 1;
b= zeros(n+1,1); b(n+1) = 1;
m = A\b; % get characteristic poly coeffs
end