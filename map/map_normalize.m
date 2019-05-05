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

for i=1:size(MAP{1},1)
    for j=1:size(MAP{1},1)
        MAP{1}(i,j)=real(MAP{1}(i,j));
        MAP{2}(i,j)=real(MAP{2}(i,j));
    end
end
MAP{1}(find(MAP{1}<0))=0;
MAP{2}(find(MAP{2}<0))=0;
for n=1:size(MAP{1},1)
    MAP{1}(n,n)=0;
    for b=1:length(MAP)
        MAP{1}(n,n)=MAP{1}(n,n)-sum(MAP{b}(n,:));
    end
end
end
