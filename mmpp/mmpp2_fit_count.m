function mmpp = mmpp2_fit_count(mu, bt1, bt2, binf, m3t2, t1, t2)
% Fits a MMPP according to [Heffes and Lucantoni, 1986].
% mu: arrival rate
% bt1: IDC at scale t1
% bt2: IDC at scale t2
% binf: IDC for t->inf
% m3t2: third central moment
% t1
% t2

if abs(binf-1) < 1e-8 && abs(binf-bt1) < 1e-8
    % degenerate case, r2 = 0
    % therefore never goes back to state 1
    % fit a marked poisson process
    %fprintf('Warning: marked Poisson process, constant IDC(t) = 1\n');
    mmpp = {-mu,mu};
    return;
end

if ~(binf>bt1 && bt1>1)
    %fprintf('Warning: no valid MMPP(2), infeasible IDC(t): IDC(%.2f) = %.3f, IDC(inf) = %.3f\n', t1, bt1, binf);
    %mmpp = {};
    mmpp = {-mu,mu};
    return;
end

% d = r1 + r2
c = (binf-1)/(binf-bt1);
% according to WolframAlpha
% the solution can be written as (ProductLog[-c e^(-c)]+c)/t1
z = -c*exp(-c);
w = fsolve(@(w) z-w*exp(w),  1, optimset('Display','none',...
                                         'MaxIter',10000,...
                                         'MaxFunEvals',10000,...
                                         'TolFun', 1.0e-12,...
                                         'TolX',1.0e-12));
x = (w + c)/t1;

k1 = mu^3 * t2^3;
k2 = 3*mu^2*(binf-1)*t2^2;
k3 = 3*mu*(binf-1)/x*t2;
k4 = 3*mu/x^2*(binf-1)*t2*exp(-x*t2);
k5 = 6*mu/x^3*(binf-1)*(1-exp(-x*t2));
g1t2 = m3t2 + 3*mu*t2*(mu*t2-1)*bt2 + mu*t2*(mu*t2-1)*(mu*t2-2);
h = (g1t2 - k1 - k2 - k3*(-mu) - k4*mu*x) / ((k3/x) + k4 -k5);
if abs(h) < 1e-4 % h=0 case
    r1 = x/2;
    r2 = x/2;
    l2 = mu - 1/2 * sqrt(2*(binf-1)*mu*x);
    l1 = mu + 1/2 * sqrt(2*(binf-1)*mu*x);
else
    y = (binf-1)*mu*x^3/(2*h^2);
    r1 = x/2 * (1 + 1/sqrt(4*y + 1));
    r2 = x - r1;
    if (r1 < r2)
        tmp = r1;
        r1 = r2;
        r2 = tmp;
    end
    w = h/(r1-r2);
    w_min = -mu/r1*(r1+r2);
    w_max = mu/r2*(r1+r2);
    if (w < w_min || w > w_max)
        %fprintf('Warning: ignoring third moment to achieve feasibility\n');
        z = (binf-1)*x^3*mu;
        u = x*z/(2*mu^2*x^2+z);
        r1 = u+(x-u)/2;
        r2 = x-r1;
        delta = sqrt(z/(2*r1*r2));
        l2 = mu - r2/x * delta;
        l1 = l2 + delta;
    else
        l2 = mu - h/(r1-r2)*(r2/(r1+r2));
        l1 = h/(r1-r2) + l2;
    end
end

mmpp = {[-(r1+l1) r1; r2 -(r2+l2)], [l1 0; 0 l2]};

end