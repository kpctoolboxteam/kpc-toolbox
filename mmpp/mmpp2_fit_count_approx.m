function [FIT] = mmpp2_fit_count_approx(a, bt1, bt2, binf, m3t2, ...
    t1, t2)
% Fits a second-order Marked MMPP.
% a: arrival rate
% bt1: IDC at scale t1
% bt2: IDC at scale t2
% binf: IDC for t->inf
% m3t2: third central moment
% t1: first time scale
% t2: second time scale

x1 = optimvar('x1'); % l1
x2 = optimvar('x2'); % l2
x3 = optimvar('x3'); % r1
x4 = optimvar('x4'); % r2

obj = fcn2optimexpr(@compute_obj, x1,x2,x3,x4, 'OutputSize', [1,1]);
prob = optimproblem('Objective',obj);

prob.Constraints.c1 = x1>=1e-6;
prob.Constraints.c2 = x2>=1e-6;
prob.Constraints.c3 = x3>=0;
prob.Constraints.c4 = x4>=0;

guess = struct;
guess.x1 = a*3/4;
guess.x2 = a*3/2;
guess.x3 = 1/3;
guess.x4 = 2/3;
[xopt, fxopt] = solve(prob, guess);

%fprintf('Unified fitting error: %f\n', fxopt);
FIT = assemble_mmap(xopt.x1, xopt.x2, xopt.x3, xopt.x4);
% elseif strcmp(method_d0d1,'gs_fmincon') == 1
%     problem = createOptimProblem('fmincon', ...
%                                  'objective', @compute_obj, ...
%                                  'x0', x0, ...
%                                  'lb', lb, ...
%                                  'ub', ub, ...
%                                  'options', options);
%     [xopt,fxopt] = run(GlobalSearch,problem);
%     %fprintf('Unified fitting error: %f\n', fxopt);
%     FIT = assemble_mmap(xopt);
% elseif strcmp(method_d0d1,'exact') == 1
%     FIT = mmpp2_fit_count(a, bt1, bt2, binf, m3t2, t1, t2);
%     if isempty(FIT)
%         % this occurs when IDC(t) < 1
%         error('Feasibility cannot be repaired.');
%     end
% else
%     error('Invalid method ''%s''\n', method);
% end

%fa = map_count_mean(FIT,1);
%fprintf('Forcing exact rate: from %f to %f\n', fa, a);
FIT = map_scale(FIT, 1/a);
fmt2 = map_count_moment(FIT,t2,1:3);

% fa = map_count_mean(FIT,1);
% fbt1 = map_count_var(FIT,t1)/map_count_mean(FIT,t1);
% fbt2 = map_count_var(FIT,t2)/map_count_mean(FIT,t2);
% fbinf = map_count_var(FIT,1e8)/map_count_mean(FIT,1e8);
% fm3t2 = fmt2(3) - 3*fmt2(2)*fmt2(1) + 2*fmt2(1)^3;
%fprintf('Rate: input = %.3f, output = %.3f\n', a, fa);
%fprintf('IDC(t1): input = %.3f, output = %.3f\n', bt1, fbt1);
%fprintf('IDC(t2): input = %.3f, output = %.3f\n', bt2, fbt2);
%fprintf('IDC(inf): input = %.3f, output = %.3f\n', binf, fbinf);
%fprintf('M3(t2): input = %.3f, output = %.3f\n', m3t2, fm3t2);

   function obj = compute_obj (l1,l2,r1,r2)
        % compute characteristics
        xa = (l1*r2 + l2*r1)/(r1 + r2);
        factor = a/xa;
        
        xbt1 = (r1*(2*l1^2*r2^2*t1*factor - 2*l2^2*r2 - 2*l1^2*r2 + 2*l2^2*r2^2*t1*factor + 4*l1*l2*r2 + 2*l1^2*r2*exp(- r1*t1*factor - r2*t1*factor) + 2*l2^2*r2*exp(- r1*t1*factor - r2*t1*factor) - 4*l1*l2*r2^2*t1*factor - 4*l1*l2*r2*exp(- r1*t1*factor - r2*t1*factor)) + r1^2*(2*r2*t1*factor*l1^2 - 4*r2*t1*factor*l1*l2 + 2*r2*t1*factor*l2^2))/(t1*factor*(r1 + r2)^3*(l1*r2 + l2*r1)) + 1;
        
        if (t1 ~= t2)
            xbt2 = (r1*(2*l1^2*r2^2*t2*factor - 2*l2^2*r2 - 2*l1^2*r2 + 2*l2^2*r2^2*t2*factor + 4*l1*l2*r2 + 2*l1^2*r2*exp(- r1*t2*factor - r2*t2*factor) + 2*l2^2*r2*exp(- r1*t2*factor - r2*t2*factor) - 4*l1*l2*r2^2*t2*factor - 4*l1*l2*r2*exp(- r1*t2*factor - r2*t2*factor)) + r1^2*(2*r2*t2*factor*l1^2 - 4*r2*t2*factor*l1*l2 + 2*r2*t2*factor*l2^2))/(t2*factor*(r1 + r2)^3*(l1*r2 + l2*r1)) + 1;
        end
        
        xbinf = ((2*r2*l1^2 - 4*r2*l1*l2 + 2*r2*l2^2)*r1^2 + (2*l1^2*r2^2 - 4*l1*l2*r2^2 + 2*l2^2*r2^2)*r1)/((r1 + r2)^3*(l1*r2 + l2*r1)) + 1;
        
        t = t2*factor;
        d = r1+r2;
        p = (l1-l2)*(r1-r2);
        xg3t = xa^3*t^3 + ...
            3*xa^2*(xbinf-1)*t^2 + ...
            3*xa*(xbinf-1)/d*(p/d-xa)*t + ...
            3*xa/d^2*(xbinf-1)*(p + xa*d)*t*exp(-t*d) - ...
            6*xa/d^3*(xbinf-1)*p*(1-exp(-t*d));
        xm3t2 = xg3t - 3*xa*t*(xa*t-1)*xbt2 - xa*t*(xa*t-1)*(xa*t-2);
        % compute objective
        obj = 0;
        obj = obj + (xa/a-1)^2;
        obj = obj + (xbt1/bt1-1)^2;
        if t1 ~= t2
            obj = obj + (xbt2/bt2-1)^2;
        end
        obj = obj + (xbinf/binf-1)^2;
        obj = obj + (xm3t2/m3t2-1)^2;
        
    end


  
    function mmap = assemble_mmap(l1,l2,r1,r2)
        % assemble
        D0 = [-(l1+r1), r1; r2, -(l2+r2)];
        D1 = [l1 0; 0 l2];
        mmap = {D0, D1};
    end

end
