function [TF] = map_issym(MAP)
% Returns 1 if the MAP or MMAP is symbolic, 0 otherwise.

TF = 0;

for i = 1:length(MAP)
    if strcmp(class(MAP{1}),'sym')
    	TF = 1;
    	return
    end
end