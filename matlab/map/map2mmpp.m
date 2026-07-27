function [Q,LAMBDA]=map2mmpp(MAP)
if (norm(MAP{2}-diag(diag(MAP{2})))>1e-10)
    warning('The MAP is not a MMPP, LAMBDA is not diagonal')
end
Q=MAP{1}+MAP{2};
LAMBDA=MAP{2};
end
