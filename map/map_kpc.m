function MAP=map_kpc(MAPa,MAPb)
% MAP=map_kpc(MAPa,MAPb)
% Kronecker product composition of MAPa and MAPb
if nargin == 1 && iscell(MAPa)
    K=length(MAPa);
    MAP = map_kpc(MAPa{1},MAPa{2});
    for k=3:K
        MAP = map_kpc(MAP,MAPa{k});
    end
else
    MAP={-kron(MAPa{1},MAPb{1}),kron(MAPa{2},MAPb{2})};
end
end