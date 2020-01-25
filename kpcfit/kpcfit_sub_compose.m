function [kpcMAP,subMAPs,error]=kpcfit_sub_compose(E1j,SCVj,E3j,G2j,varargin)
% map_kpcfit_compose make a large MAP out of several small MAP(2) or
% MA(1)s.
%
%  MAP=map_kpcfit_compose(E1j,SCVj,E3j,G2j) gives a big MAP composed by
%  many small MAP(2) or MAP(1)s with individual E1, SCV, E3, G2.
%
%   Input:
%       E1j     = a vector of the first moment for every map
%       SCVj    = a vector of squared coefficient of variation for every
%       map
%       E3j     = a vector of the third moment for every map
%       G2j     = a vector of gamma2 for every map
%
%   Output:
%       MAP     = A map composed by J maps, each map(2) has the specified
%       first moment, SCV, third moment, gamma2 in the vectors of E1j,
%       SCVj, E3j and so on.
%
%   Comment:
%       The first map is a mmpp2 while the other J-1 maps are HMAPs
%
J=length(G2j);

verbose=0;
if nargin >= 5
    verbose = varargin{1};
end

% compute the composed process
error=0;
kpcMAP = mmpp2_fit3(E1j(1),(1+SCVj(1))*E1j(1)^2,E3j(1),G2j(1));
if any(imag(kpcMAP{1})>10^-4, 'all') || any(imag(kpcMAP{2})>10^-4, 'all') || map_isfeasible(kpcMAP) == 0
    if SCVj(1)<0.5
        kpcMAP = map_erlang(E1j(1),2);
        if verbose > 0
            fprintf(1,'MAP 1 is erlang-2\n');
        end
    else
        if verbose > 0
            fprintf(1,'MAP 1 has presumably infeasible E3: %f\n', E3j(1));
        end
        [kpcMAP,fitErr] = map2_fit(E1j(1),(1+SCVj(1))*E1j(1)^2,-1,G2j(1));
        if isempty(kpcMAP) || fitErr ~= 0
           if verbose > 0
            fprintf(1,'MAP 1 has infeasible G2: %f (err=%d)\n', G2j(1), fitErr);
           end
           [kpcMAP,fitErr] = map2_fit(E1j(1),(1+SCVj(1))*E1j(1)^2,-1,0);
           
           if fitErr ~= 0 || isempty(kpcMAP)
                %fprintf("map2_fit returned: %d\n", fitErr);
                error=1;
                subMAPs=cell(1,0);
                return
           end
        end
    end
end

subMAPs{1}=kpcMAP;

for j=2:J
    MAPj=map_feasblock(E1j(j),(1+SCVj(j))*E1j(j)^2,E3j(j),G2j(j));
    try
        feasible = map_isfeasible(MAPj);
    catch
        if verbose > 0
            fprintf(1, "Could not check feasibility of MAP-%d", j);
        end
        feasible = 0;
    end
    if feasible == 0
        if SCVj(j)<1
            if verbose > 0
                fprintf(1,'MAP %d has low variability\n',j);
            end
            MAPj = map2_exponential(E1j(j));
        else
            if verbose > 0
                fprintf(1,'MAP %d has presumably infeasible E3\n',j);
            end
            [MAPj,fitErr] = map2_fit(E1j(j),(1+SCVj(j))*E1j(j)^2,-1,G2j(j));
            if fitErr ~= 0 || isempty(MAPj)
                [MAPj,fitErr] = map2_fit(E1j(j),(1+SCVj(j))*E1j(j)^2,-1,0);
                
                if verbose > 0
                    fprintf(1,'MAP %d has infeasible G2: %f)\n',j, G2j(j));
                end
                
                if fitErr ~= 0 || isempty(MAPj)
                    error=5;
                    return
                end
            end            
            lambda=eig(-inv(MAPj{1}));
            D0 = diag(-1./(lambda));
            P = map_embedded(MAPj);
            p = min(eig(P));
            D1 = -D0*[p,1-p;p,1-p];
            MAPj=map_normalize({D0,D1});
        end
    end
    if isempty(MAPj)
        if verbose > 0
            fprintf(1, 'Replacing MAP %d with an exponential\n', j);
        end
        subMAPs{j}=map_exponential(E1j(j));
    else
        subMAPs{j}=MAPj;
    end
    
    kpcMAP= map_kpc(kpcMAP,MAPj);
end
for j=1:J
    if map_isfeasible(subMAPs{j}) == 0
        fprintf(1,'MAP %d is infeasible\n',j);
        error=10;
    end
end


kpcMAP=map_normalize(kpcMAP);

% for j=1:J
%     fprintf("MAP %d: ", j);
%     fprintf("E1: %f (%f)\n", map_mean(subMAPs{j}), E1j(j));
%     fprintf("SCV: %f (%f)\n", map_scv(subMAPs{j}), SCVj(j));
%     fprintf("G2: %f (%f)\n", map_gamma(subMAPs{j}), G2j(j));
% end
end
