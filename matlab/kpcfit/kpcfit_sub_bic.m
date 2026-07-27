function Nstar =  map_kpcfit_bic(SA,orders,varargin)
% function order = map_kpcfit_bic(SA,orders,...) returns the possible order for
% a trace for the purpose of KPCfitting
%
%   order = map_kpcfit_bic(SA,orders,...) returns the order of KPCmap as a power of 2,
%   SA is a sequence of autocorrelations of a trace, ordered from lag-1 to lag-n, and
%   orders is a range of orders that will be chosen from.


nlags = length(SA);
nlagsend = find(SA < 1e-6); % nlagsend = point of first decay below zero
ordermax = max(orders);
if isempty(nlagsend)
    nlagsend = nlags;
else
    nlagsend = nlagsend(1)-1 + ordermax;
end

NLAGSMAX = 10000;

if (nlagsend > NLAGSMAX)
    SAlags = unique(round(logspace(log10(1),log10(nlagsend-ordermax),NLAGSMAX)));
    len = length(SAlags);
    temp = NLAGSMAX - length(SAlags);
    step = (nlagsend-ordermax)^(1/NLAGSMAX);
    i = 1;
    while (((SAlags(i) - round(step^i)) < temp) && (i<=len))
        i = i + 1;
        
    end
    SAlags = [1:SAlags(i),SAlags(i+1:end)];
else
    % if nlagsend-ordermax*2 < 0, suggest something to the user
    if nlagsend > (ordermax)
        SAlags = 1:(nlagsend-ordermax);
    else
        ordermax = nlagsend-2;
        SAlags = 1:(nlagsend-ordermax);
        orders = orders(find(orders<=ordermax));
    end
end
%SAlags = linspace(ordermax+1,ordermax+1000,500);

n_samples = length(SAlags);
SAlags = round(SAlags);
X = [];
X = SA(SAlags);
rarray = [];
parray = [];
rsarray = [];
adjRarray = [];
cp = [];
AIC = [];
mean_error = [];
predict_error = [];
SBC = [];
sse_sacf = [];
bpre = zeros(2,1);
SAlags_y = SAlags;
%SAlags_y = SAlags + ordermax;
Y = SA(SAlags_y);
j = ordermax;
order = j;
X = [];
for i = 1:order
    lags_x = SAlags_y + i;
    X(:,end+1) = SA(lags_x);
end

Y_var = var(Y);
%disp(size(Y));
%disp(size(X));
%X(:,end+1) = ones(size(X,1),1);
warning off
if length(Y) == 1
    %[b,bint,r,rint,stats] = regress([Y;Y],[X;X])
    MSEtotal = 1e6;
else
    [b,bint,r,rint,stats] = regress(Y,X);
    MSEtotal = sum(r.^2)/(n_samples-1);
end
warning on


for j = 1:length(orders)
    order = orders(j);
    X = [];
    for i = 1:order
        lags_x = SAlags_y + i;
        X(:,end+1) = SA(lags_x);
    end
    
    Y_var = var(Y);
    %X(:,end+1) = ones(size(X,1),1);
    warning off all;
    [b,bint,r,rint,stats] = regress(Y,X);
    %predict_error(end+1) = norm(Ypredict - Xpredict(:,1:order)*b);
    AIC(end+1) = n_samples*log(sum(r.^2)) - n_samples*log(n_samples) + 2*order;
    SBC(end+1) = n_samples*log(sum(r.^2)) - n_samples*log(n_samples) + log(n_samples)*order;
    cp(end+1) = sum(r.^2)/MSEtotal - (n_samples-2*order);
    mean_error(end+1) = mean(r);
    sse_sacf(end+1) = sum(r.^2)/sum(Y.^2);
    %disp(bint);
    adjR = 1 - sum(r.^2)/(n_samples-order)/Y_var;
    adjRarray(end+1) = adjR;
    bpre = b;
    r = norm(r,2);
    rarray(end+1) = r;
    parray(end+1) = stats(4);
    rsarray(end+1) =  stats(1);
    
    
    %disp(stats(3));
end
xaxis = orders;
%h0 = plot(xaxis,mean_error);
%h1 = plot(xaxis,rarray);
%h2 = plot(xaxis,rsarray);
%h3 = plot(xaxis,adjRarray);
%h3 = plot(xaxis,parray);
%h5 = plot(xaxis,cp);
%h6 = plot(xaxis,AIC);
%h7 = plot(xaxis,predict_error);
if ~isempty(varargin)
    
    h8 = plot(xaxis,SBC,varargin{1});
end
%h9 = plot(xaxis,sse_sacf);
%hold on;
%h2 = plot(xaxis,parray,'*');

%h = gca;
%legend(h,'R-norm','p-value');


[v,ind] = sort(SBC,'ascend');

Nstar =  log2(orders(ind(1)));

order = 2^Nstar;
X = [];
for i = 1:order
    lags_x = SAlags_y + i;
    X(:,end+1) = SA(lags_x);
end

warning off all;

[b,bint,r,rint,stats] = regress(Y,X);
rssrho = sum(r.^2)/sum(Y.^2);

if (rssrho > 0.05)
    disp('Warning! The residual sum of square indicates that the trace may not be well characterized by a MAP with less than 64 states!');
    
    disp('##################################### Using another model for fitting is suggested! ########################################');
    
end



% a small example for the order 16 case


end
