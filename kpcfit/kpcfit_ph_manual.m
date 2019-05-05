function PH = kpcfit_ph_manual(E,varargin)
%  ** beta version **
%
%  PH = kpcfit_ph_manual(E,'option1','val1','option2','val2',...)
%
% DESCRIPTION
% Fit phase-type (PH) process using Kronecker Product Composition (KPC)
% method
%
% INPUT
%
% E            - vector of moments of consecutive order to be fitted
%
% EXAMPLE
%
% S % trace
% E = [];
% for k=1:6 E(k)=mean(S.^k); end
% PH = kpcfit_ph_manual(E);
%
% OPTION LIST
% 'WantHyper'    - boolean, if true returns hyper-exponential distribution
% 'NumStates'  - integer, number of states of the PH distribution - **must be multiple of 2**, default: 4
% 'NumPHs'     - integer, equals log2(number of states of the PH distribution), default: 2 (i.e., PH(4))
% 'NumMoms'    - integer, moments considered in the moment matching problem, default: 6
% 'MaxIterMM'  - integer, maximum number of iterations for a single moment-matching fitting run

%% options
OptionNames = [
    'WantHyper ';
    'NumStates ';
    'NumPHs    ';
    'NumMoms   ';
    'MaxIterMM ';
    ];

OptionTypes = [
    'numeric';
    'numeric';
    'numeric';
    'numeric';
    'numeric'];

OptionValues = [];
for i = 1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end

% Default settings
options.NumStates = 4;
options.NumPHs = 2;
options.NumMoms = 6;
options.MaxIterMM = 100;

% Parse optional parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);
if mod(options.NumStates,2) ~= 0
    error('kpcfit_ph: number of states must be multiple of 2');
end
%
TOL=1e-10;
EPSTOL=100*TOL;
warning off
optimoptions = optimset('Algorithm','active-set', ...
    'LargeScale','off', ...
    'MaxIter',options.MaxIterMM, ...
    'MaxFunEvals',1e10,  ...
    'MaxSQPIter',500,  ...
    'TolCon',TOL,  ...
    'Display','off', ...
    'OutputFcn',@kpcfit_ph_manual_outfun);

stagnval = 0;
stagniter = 0;
x0=zeros(3*options.NumPHs,1); % each PH brings 3 degrees of freedomg
lgkx = x0; % last good known x
n=length(E); % number of moments
nexact=2;

for j=1:options.NumPHs
    if j==1
        PH=map_rand(2);
    else
        PH=map_feasblock(rand,1000*rand,-1,0); % hyperexponential
    end
    for i=1:3
        x0((i-1)*options.NumPHs+j)=map_moment(PH,i);
    end
end

fmincon(@objfun,x0,[],[],[],[],x0*0+EPSTOL,[],@nnlcon,optimoptions);
xtopar(lgkx); % this recomputes PH
PH=map_normalize(map_scale(PH,E(1)));

if map_isfeasible(PH)<10
    warning('PH distribution may not be valid, check for negative rates in D1 or in off-diagonal of D0')
end
PH=map_scale(PH,E(1));


    function [Eest]=xtopar(x)
        x=abs(x)+EPSTOL;
        
        % extract moments
        E1j=x(1:options.NumPHs)';
        SCVj=x((options.NumPHs+1):2*options.NumPHs)';
        E2j=(1+SCVj).*(E1j.^2);
        E3j=x((2*options.NumPHs+1):3*options.NumPHs);
        
        % fit ph
        if options.WantHyper
            PH = map_feasblock(E1j(1),E2j(1),E3j(1),rand);
        else
            PH = mmpp2_fit3(E1j(1),E2j(1),E3j(1),rand);
        end
        PH = map_scale(PH,E1j(1));

        for j=2:options.NumPHs-1
            [PHk]=map_feasblock(E1j(1),E2j(j),E3j(j),rand);
            PH=map_kpc(PH,PHk);
        end
        E1j(options.NumPHs)=E(1)/map_moment(PH,1);
        E2j(options.NumPHs)=E(2)*2/map_moment(PH,2);
        if length(E)>2
            E3j(options.NumPHs)=E(3)*6/map_moment(PH,3);
        end
        x(2)=E1j(options.NumPHs);
        x(4)=E2j(options.NumPHs);
        x(6)=E3j(options.NumPHs);
%        map_scv(PH)
%        map_isfeasible(PH)
        [PH]=map_kpc(PH,map_feasblock(E1j(options.NumPHs),E2j(options.NumPHs),E3j(options.NumPHs),0));
%        F=map_feasblock(E1j(options.NumPHs),E2j(options.NumPHs),E3j(options.NumPHs),0);
%         F{1}
%         F{2}
%         %map_scv(PH)
%         map_isfeasible(PH)
%         %PH{1}        
%        pause
        
        PH{1}=diag(diag(PH{1}));
        for k=1:length(E)
            Eest(k)=map_moment(PH,k);
        end
    end

    function [c,ceq]=nnlcon(x)
        [Eest]=xtopar(x);
        c=[];
        ceq=[];
        for k=3:nexact
            ceq(end+1)=norm((E(k))-(Eest(k)),1);
        end
    end

    function f=objfun(x)
        [Eest]=xtopar(x);
        f=0;
        for k=(nexact+1):length(E)
            f=f+norm(log10(E(k))-log10(Eest(k)),1);
        end
    end

end

function [MAP]=map_rand(K)
if nargin<1
    K=2;
end
D1=rand(K,K);
D0=rand(K,K);

MAP=cell(1,2);
MAP{1}=D0;
MAP{2}=D1;
MAP=map_normalize(MAP);
end