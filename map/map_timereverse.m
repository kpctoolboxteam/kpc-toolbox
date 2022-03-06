function MAPr=map_reverse(MAP)
piq=map_pi(MAP);
D=diag(piq);
MAPr={inv(D)*MAP{1}'*D,inv(D)*MAP{2}'*D};
end
