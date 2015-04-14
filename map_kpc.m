function MAP=map_kpc(MAPa,MAPb)
% MAP=map_kpc(MAPa,MAPb)
% Kronecker product composition of MAPa and MAPb
    MAP={-kron(MAPa{1},MAPb{1}),kron(MAPa{2},MAPb{2})};
end