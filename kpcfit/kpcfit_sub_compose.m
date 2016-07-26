function [MAP,subMAPs]=kpcfit_sub_compose(E1j,SCVj,E3j,G2j)
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
% compute the composed process

MAP = mmpp2_fit3(E1j(1),(1+SCVj(1))*E1j(1)^2,E3j(1),G2j(1));
if map_isfeasible(MAP) == 0
    if SCVj(1)<0.5
        MAP = map_erlang(E1j(1),2);
%        fprintf(1,'MAP 1 is erlang-2\n');
    else
        MAP = map2_fit(E1j(1),(1+SCVj(1))*E1j(1)^2,-1,G2j(1));
%        fprintf(1,'MAP 1 has presumably infeasible E3\n');
        if isempty(MAP)
            MAP = map2_fit(E1j(1),(1+SCVj(1))*E1j(1)^2,-1,0);
%           fprintf(1,'MAP 1 has infeasible G2\n');
        end
    end
end
subMAPs{1}=MAP;
for j=2:J
    MAPj=map_feasblock(E1j(j),(1+SCVj(j))*E1j(j)^2,E3j(j),G2j(j));
    if map_isfeasible(MAPj) == 0
        if SCVj(j)<1
            MAPj = map2_exponential(E1j(j));
%            fprintf(1,'MAP has low variability\n',j);
        else
            MAPj = map2_fit(E1j(j),(1+SCVj(j))*E1j(j)^2,-1,G2j(j));
%            fprintf(1,'MAP %d has presumably infeasible E3\n',j);
            if isempty(MAPj)
                MAPj = map2_fit(E1j(j),(1+SCVj(j))*E1j(j)^2,-1,0);
%                fprintf(1,'MAP %d has infeasible G2\n',j);
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
        subMAPs{j}=map_exponential(E1j(j));
    else
        subMAPs{j}=MAPj;
    end
    MAP= map_kpc(MAP,MAPj);
end
for j=1:J
    if map_isfeasible(subMAPs{j})
%        fprintf(1,'MAP %d is feasible\n',j);
    else
%        fprintf(1,'MAP %d is infeasible\n',j);
    end
end
MAP=map_normalize(MAP);
end
