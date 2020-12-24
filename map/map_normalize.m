function MAP=map_normalize(MAP)
% MAPOUT=map_normalize(MAPIN) - Try to make a MAP feasible
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

D0=MAP{1};
D1=MAP{2};
for i=1:size(D1,1)
    for j=1:size(D1,1)
        D0(i,j)=real(D0(i,j));
        D1(i,j)=real(D1(i,j));
    end
end
D0(find(D0<0))=0;
D1(find(D1<0))=0;
for n=1:size(D0,1)
    D0(n,n)=0;
    for b=1:length(MAP)
        D0(n,n)=D0(n,n)-sum(MAP{b}(n,:));
    end
end
MAP = {D0,D1};
end
