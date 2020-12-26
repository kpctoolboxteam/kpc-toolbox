function MAP = map_normalize(MAP)
% MAPOUT = map_normalize(MAPIN) - Try to make a MAP feasible
%
%  Input:
%  MAPIN: a MAP in the form of {D0,D1}
%
%  Output:
%  MAPOUT: MAPIN with D0 normalized to make D0+D1 an infinitesimal
%  generator, with all negative entries set to zero, and with complex
%  values set equal to their norm
%
%  Examples:
%  - map_normalize({[0,0;0,0],[1,2;3,4]}) produces a valid MAP
%

D0=real(MAP{1});
D1=real(MAP{2});
D0 = D0 - diag(diag(D0)); % remove diagonal
D0(D0<0)=0; % remove invalid entries
D0 = D0 - diag(sum(D0));
MAP{1} = D0;
for b=2:length(MAP) %% D1..Dn
    Db = MAP{b};
    Db(Db<0)=0;
    D0 = D0 - diag(sum(Db,2));
    MAP{b} = Db;
end
MAP{1} = D0;
end
