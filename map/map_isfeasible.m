function ISFEAS=map_isfeasible(MAP, TOL)
% ISFEAS=map_isfeasible(MAP) - Evaluate feasibility of a MAP process
%
%  Input:
%  MAP: a MAP in the form of {D0,D1}
%
%  Output:
%  ISFEAS: boolean 1=feasible, 0=infeasible. Numerical tolerance is based
%  on the standard toolbox value in map_feastol.m
%
%  Examples:
%  - map_isfeasible({[0,0;0,0],[1,2;3,4]}) is an infeasible MAP
%

ISFEAS=0;
if nargin==2
    ISFEAS=map_checkfeasible(MAP,TOL);
elseif nargin==1
    TOLMAGNITUDE=15;
    for k=TOLMAGNITUDE:-1:1
        check=map_checkfeasible(MAP,10^-k);
        if check
            ISFEAS=k;
            ISFEAS=ISFEAS>map_feastol;
            return
        end
    end
    if nargout>0
        ISFEAS=0;
    end
end
end

function isfeas=map_checkfeasible(MAP,TOL)
w=warning;
warning off;
n=size(MAP{1});
D0=MAP{1};
D1=MAP{2};
if sum(isnan(D0(:)))>0
    warning(w);
    isfeas=0;
    return
end

if any(isinf(D0), 'all') || any(isinf(D1), 'all') || any(isnan(D0), 'all') || any(isnan(D1), 'all')
    isfeas=0;
    return
end
if any(imag(D0)>10^-4, 'all') || any(imag(D1)>10^-4, 'all') 
    isfeas=0;
    return
end


P=inv(-D0)*D1;
Q=D0+D1;

% set very small values to zero
D0(find(abs(D0)<TOL))=0;
D1(find(abs(D1)<TOL))=0;
P(find(abs(P)<TOL))=0;
Q(find(abs(Q)<TOL))=0;

% assume initially the MAP is feasible
isfeas=1;
for i=1:n
    for j=1:n
        %% validate signs
        if i~=j && D0(i,j)<0
            isfeas=0;
        end
        if i==j && D0(i,j)>0
            isfeas=0;
        end
        if D1(i,j)<0
            isfeas=0;
        end
        if i~=j && Q(i,j)<0
            isfeas=0;
        end
        if i==j && Q(i,j)>0
            isfeas=0;
        end
        if P(i,j)<0
            isfeas=0;
        end
    end
    %% validate stochasticity
    if abs(sum(P(i,:)))<1-n*TOL
        isfeas=0;
    end
    if abs(sum(P(i,:)))>1+n*TOL
        isfeas=0;
    end
    if abs(sum(Q(i,:)))<0-n*TOL
        isfeas=0;
    end
    if abs(sum(Q(i,:)))>0+n*TOL
        isfeas=0;
    end
end
if isfeas==0
    warning(w);
    return;
end
%% validate irreducibility
if n<map_largemap % compute eigenvalues only if the MAP is not too large
    eigQ=eig(Q);
    if length(find(eigQ>0-TOL))>1
        isfeas=0;
    end
    eigP=eig(P);
    if length(find(eigP>1-TOL))>1
        isfeas=0;
    end
end
end
