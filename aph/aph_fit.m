function APH = aph_fit(e1,e2,e3,nmax)
% APH = aph_fit(e1,e2,e3,nmax) - moment matching of APH(n) distribution
% nmax = max order (default: 10)
%
% Implementation of: A.Bobbio, A.Horvath, M.Telek, Matching three moments
% with minimal acyclic phase type distributions, Stochastic Models
% 21:303-326, 2005.

if nargin < 4
    nmax = 10;
end
if isinf(e2) || isinf(e3)
    APH = map_exponential(e1);
    return
end
n2 = e2/e1^2;
n3 = e3/e1/e2;

% find order
n2_feas = false;
n3_ubfeas = false;
n3_lbfeas = false;
n = 1;
un = 0;
while (n2_feas == false || n3_lbfeas == false || n3_ubfeas == false ) && n < nmax
    n = n + 1;
    pn = ((n+1)*(n2-2)/(3*n2*(n-1)))*(-2*sqrt(n+1)/sqrt(4*(n+1)-3*n*n2) - 1);
    an = (n2 - 2) / (pn*(1-n2) + sqrt(pn^2+pn*n*(n2-2)/(n-1)));
    ln = ((3+an)*(n-1)+2*an)/((n-1)*(1+an*pn)) - (2*an*(n+1))/(2*(n-1)+an*pn*(n*an+2*n-2));
    un_1 = un;
    un = (1/(n^2*n2))*(2*(n-2)*(n*n2-n-1)*sqrt(1+n*(n2-2)/(n-1))+(n+2)*(3*n*n2-2*n-2));
    if n2 >= (n+1)/n && n2 <= (n+4)/(n+1)
        n2_feas = true;
        if n3 >= ln
            n3_lbfeas = true;
        end
    elseif n2 >= (n+4)/(n+1)
        n2_feas = true;
        if n3 >= n2*(n+1)/n
            n3_lbfeas = true;
        end
    end
    if n2 >= (n+1)/n && n2 <= n/(n-1)
        n2_feas = true;
        if n3 <= un
            n3_ubfeas = true;
        end
    elseif n2 >= n/(n-1)
        n2_feas = true;
        if n3 < Inf
            n3_ubfeas = true;
        end
    end
end
%keyboard
if (n2_feas == false || n3_lbfeas == false || n3_ubfeas == false ) || (n == nmax)
    warning('cannot match moment set exactly');
    n2 = (n+1)/n;
    n3 =  2*n2-1;
    %n3 = 2*n2-1;
    %n3 = un
end

% fitting
if n2 <= n/(n-1) || n3 <= 2*n2-1 % case 1 of 2
    b = 2*(4-n*(3*n2-4))/(n2*(4+n-n*n3)+sqrt(n*n2)*sqrt(12*n2^2*(n+1)+16*n3*(n+1)+n2*(n*(n3-15)*(n3+1)-8*(n3+3))));
    a = (b*n2-2)*(n-1)*b/((b-1)*n);
    p = (b-1)/a;
    lambda = 1;
    mu = lambda*(n-1)/a;
    alpha = zeros(1,n); alpha(1)=p; alpha(end)=1-p;
    T = diag(-mu*ones(1,n))+diag(mu*ones(1,n-1),1); T(end,end)=-lambda;
    APH = map_scale(map_normalize({T,-T*ones(n,1)*alpha}),e1);
elseif n2 > n/(n-1) && n3 > un_1 % case 2 of 2
    K1 = n-1;
    K2 = n-2;
    K3 = 3*n2-2*n3;
    K4 = n3-3;
    K5 = n-n2;
    K6 = 1+n2-n3;
    K7 = n+n2-n*n2;
    K8 = 3+3*n2^2+n3-3*n2*n3;
    K9 = 108*K1^2*(4*K2^2*K3*n^2*n2+K1^2*K2*K4^2*n*n2^2+4*K1*K5*(K5^2-3*K2*K6*n*n2)+sqrt(-16*K1^2*K7^6+(4*K1*K5^3+K1^2*K2*K4^2*n*n2^2+4*K2*n*n2*(K4*n^2-3*K6*n2+K8*n))^2));
    K10 = K4^2/(4*K3^2) - K5/(K1*K3*n2);
    K11 = 2^(1/3)*(3*K5^2+K2*(K3+2*K4)*n*n2)/(K3*K9^(1/3)*n2);
    K12 = K9^(1/3) / (3*2^(7/3)*K1^2*K3*n2);
    K13 = sqrt(K10 + K11 + K12);
    K14 = (6*K1*K3*K4*K5+4*K2*K3^2*n-K1^2*K4^3*n2) / (4*K1^2*K3^3*K13*n2);
    K15 = -K4/(2*K3);
    K16 = sqrt(2*K10 - K11 -K12 -K14);
    K17 = sqrt(2*K10 - K11 -K12 +K14);
    K18 = 36*K5^3 + 36*K2*K4*K5*n*n2 + 9*K1*K2*K4^2*n*n2^2 - sqrt(81*(4*K5^3+4*K2*K4*K5*n*n2+K1*K2*K4^2*n*n2^2)^2-48*(3*K5^2+2*K2*K4*n*n2)^3);
    K19 = -K5/(K1*K4*n2) -2^(2/3)*(3*K5^2+2*K2*K4*n*n2)/(3^(1/3)*K1*K4*n2*K18^(1/3)) - K18^(1/3)/(6^(2/3)*K1*K4*n2);
    K20 = 6*K1*K3*K4*K5 + 4*K2*K3^2*n - K1^2*K4^3*n2;
    K21 = K11 + K12 + K5/(2*n*K1*K3);
    K22 = sqrt(3*K4^2/(4*K3^2) - 3*K5/(K1*K3*n2) + sqrt(4*K21^2 - n*K2/(n2*K1^2*K3)));
    if n3 > un_1 && n3 < 3*n2/2
        f = K13+K15-K17;
    elseif n3 == 2*n2/2
        f = K19;
    elseif n3 > 3*n2/2 && K20 > 0
        f = -K13+K15+K16;
    elseif K20 == 0
        f = K15+K22;
    elseif K20<0
        f = K13+K15+K17;    
    end    
    
    a = 2*(f-1)*(n-1)/((n-1)*(n2*f^2-2*f+2)-n);
    p = (f-1)*a;
    lambda = 1;
    mu = lambda*(n-1)/a;
    alpha = zeros(1,n); alpha(1)=p; alpha(2)=1-p;
    T = diag(-mu*ones(1,n))+diag(mu*ones(1,n-1),1); T(1,1)=-lambda; T(1,2)=lambda;
    APH = map_scale(map_normalize({T,-T*ones(n,1)*alpha}),e1);
else
    warning('moment set cannot be matched with an APH distribution');
end

end