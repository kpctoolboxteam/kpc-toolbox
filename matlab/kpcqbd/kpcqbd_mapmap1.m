function [XN,QN,UN,pqueue] = kpcqbd_mapmap1(MAPa, MAPs)
% [XN,QN,UN,PQUEUE] = KPCQBD_MAPMAP1(MAPA,MAPS) - exact analysis of a
% MAP/MAP/1 queue via the matrix-geometric method.
%
% MAPA    - arrival process {D0,D1}
% MAPS    - service process {D0,D1}
% XN      - throughput (arrival rate)
% QN      - mean queue length (number in system)
% UN      - utilization
% PQUEUE  - queue-length distribution, one row per level, columns indexed
%           by the (arrival x service) phase; row n+1 is level n

na = length(MAPa{1});
ns = length(MAPs{1});
m  = na*ns;

A2 = kron(MAPa{2}, eye(ns));
A1 = krons(MAPa{1}, MAPs{1});
A0 = kron(eye(na), MAPs{2});
L0 = kron(MAPa{1}, eye(ns));

lambda = map_lambda(MAPa);
mu     = map_lambda(MAPs);

G = qbd_gmatrix(A0, A1, A2);
U = A1 + A2*G;
R = A2 / (-U);

M  = [L0, A2; A0, A1 + R*A0];
[~,~,V] = svd(full(M'));
xi = V(:,end)';
pi0 = xi(1:m);
pi1 = xi(m+1:end);

tot = sum(pi0) + sum(pi1 * ((eye(m)-R)\ones(m,1)));
pi0 = pi0/tot;
pi1 = pi1/tot;
if sum(pi0) < 0
    pi0 = -pi0; pi1 = -pi1;
end

maxLevels = 20000;
pqueue = zeros(maxLevels, m);
pqueue(1,:) = pi0;
row = pi1;
captured = sum(pi0);
n = 1;
while n < maxLevels
    pqueue(n+1,:) = row;
    captured = captured + sum(row);
    if captured > 1-1e-12 && n >= 2
        n = n+1; break
    end
    row = row * R;
    n = n+1;
end
pqueue = pqueue(1:n,:);

XN = lambda;
UN = 1 - sum(pqueue(1,:));
QN = (0:size(pqueue,1)-1) * sum(pqueue,2);
end
