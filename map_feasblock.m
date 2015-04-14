function MAP=map_feasblock(E1,E2,E3,G2,OPT)
if E2==2*E1^2 %exponential
    MAP=map_scale({[-1,0;0,-1],[0.5,0.5;0.5,0.5]},E1);
    return
end

% map_feasblock(E1,E2,E3,G2,OPT) - fits the most similar
if nargin==5
    if strcmpi(OPT,'scv')
        E2=(1+E2)*E1^2;
    end
else
    OPT='';
end

if E2<=2*E1^2
    warning('E2 failure (SCV=1), setting SCV=1.001');
    E2=(1+kpcfit_tol)*E1^2;
end
if E3<=(3/2)*E2^2/E1
    warning('E3 failure, setting E3=(3/2+1e-6)*E2^2/E1');
    E3=(3/2+kpcfit_tol)*E2^2/E1;
end
%[E1,E2,E3,G2]
MAP=map_block(E1,E2,E3,G2);
end